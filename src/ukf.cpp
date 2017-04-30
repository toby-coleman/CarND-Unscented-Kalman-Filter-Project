#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  // Check for initialisation
  if (!is_initialized_) {
    std::cout << "Starting initialisation" << std::endl;
    float x, y;

    // Sigma point weights
    weights_ = VectorXd(2 * n_aug_ + 1);
    weights_.fill(0.5 / (lambda_ + n_aug_));
    weights_[0] = lambda_ / (lambda_ + n_aug_);
    // State vector and covariance matrix
    x_.fill(0);
    P_ = MatrixXd::Identity(5, 5);
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Initialise with radar measurement [rho, phi, rho_dot]
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];
      x = rho * cos(phi);
      y = rho * sin(phi);
    } else {
      // Initialise with lidar measurement [x, y]
      x = meas_package.raw_measurements_[0];
      y = meas_package.raw_measurements_[1];
      // Use smaller initial covariance on x- and y-
      P_(0, 0) = 0.1;
      P_(1, 1) = 0.1;
    }
    // Don't initialise at x=0, y=0
    if ((fabs(x) < 0.0001) & (fabs(y) < 0.0001)) {
      x = 0.0001;
    }
    x_(0) = x;
    x_(1) = y;
    // Use larger initial covariance on velocity
    P_(2, 2) = 10;

    std::cout << "Initialisation complete" << std::endl;
    std::cout << "x_initial = " << x_ << std::endl;
    std::cout << "P_initial = " << P_ << std::endl;

    // Store timestamp
    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;
    return;
  }

  // Process measurement according to type
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    if (use_radar_) {
      // Predict step
      Prediction((meas_package.timestamp_ - time_us_) / 1000000.0);
      // Store timestamp
      time_us_ = meas_package.timestamp_;
      // Update
      UpdateRadar(meas_package);
    }
  } else {
    if (use_laser_) {
      // Predict step
      Prediction((meas_package.timestamp_ - time_us_) / 1000000.0);
      // Store timestamp
      time_us_ = meas_package.timestamp_;
      // Update
      UpdateLidar(meas_package);
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  // Create augmented mean vector, state covariance and sigma point matrix
  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7, 7);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug.fill(0);
  x_aug.head(5) = x_;
  P_aug.fill(0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_ * std_a_; // Q
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  // Compute sqrt of augmented P
  MatrixXd A = P_aug.llt().matrixL();
  //create augmented sigma points
  MatrixXd part = sqrt(lambda_ + n_aug_) * A;
  Xsig_aug.col(0) = x_aug;

  for (int c = 0; c < A.cols(); c++) {
      Xsig_aug.col(1 + c) = x_aug + part.col(c);
      Xsig_aug.col(1 + A.cols() + c) = x_aug - part.col(c);
  }

  // Predict sigma points
  for (int i = 0; i < Xsig_pred_.cols(); i++)
  {
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // Predicted state values
    double px_p, py_p;

    // Avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
        py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else {
        px_p = p_x + v * delta_t * cos(yaw);
        py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // Add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    // Write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  // Predict state
  x_.fill(0);
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
      x_ = x_ + weights_[i] * Xsig_pred_.col(i);
  }
  // Predict state covariance matrix
  P_.fill(0);
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
      P_ = P_ + weights_[i] * (Xsig_pred_.col(i) - x_) *
            (Xsig_pred_.col(i) - x_).transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  // Lidar measurements comprise [x, y]

  // Sigma points in measurement space
  MatrixXd Zsig = MatrixXd(2, 2 * n_aug_ + 1);
  // Mean predicted measurement
  VectorXd z_pred = VectorXd(2);
  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(2, 2);
  // Cross cross correlation matrix Tc
  MatrixXd Tc = MatrixXd(n_x_, 2);

  // Retrieve measurement, z
  VectorXd z = meas_package.raw_measurements_;

  // Transform sigma points into measurement space
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
      double p_x = Xsig_pred_(0, i);
      double p_y = Xsig_pred_(1, i);

      Zsig.col(i) << p_x, p_y;
  }

  // Calculate mean predicted measurement
  z_pred.fill(0);
  for (int i = 0; i < Zsig.cols(); i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // Calculate measurement covariance matrix S
  MatrixXd R = MatrixXd(2, 2);
  R.fill(0);
  R(0, 0) = std_laspx_ * std_laspx_;
  R(1, 1) = std_laspy_ * std_laspy_;

  S.fill(0);
  for (int i = 0; i < Zsig.cols(); i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  S = S + R;

  // Calculate cross correlation matrix
  Tc.fill(0);
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
      Tc = Tc + weights_(i) * (Xsig_pred_.col(i) - x_) *
            (Zsig.col(i) - z_pred).transpose();
  }

  // Calculate residual
  VectorXd z_diff = z - z_pred;

  // Calculate Kalman gain K;
  MatrixXd Sinv = S.inverse();
  MatrixXd K = Tc * Sinv;

  // Update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  // Calculate NIS
  NIS_laser_ = z_diff.transpose() * Sinv * z_diff;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  // Radar measurements comprise [rho, phi, rho_dot]

  // Sigma points in measurement space
  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);
  // Mean predicted measurement
  VectorXd z_pred = VectorXd(3);
  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(3, 3);
  // Cross cross correlation matrix Tc
  MatrixXd Tc = MatrixXd(n_x_, 3);

  // Retrieve measurement, z
  VectorXd z = meas_package.raw_measurements_;

  // Transform sigma points into measurement space
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
      double p_x = Xsig_pred_(0, i);
      double p_y = Xsig_pred_(1, i);
      double v = Xsig_pred_(2, i);
      double psi = Xsig_pred_(3, i);

      double rho = sqrt(p_x * p_x + p_y * p_y);
      double phi = atan2(p_y, p_x);
      double rho_dot = (p_x * cos(psi) * v + p_y * sin (psi) * v) / rho;

      Zsig.col(i) << rho, phi, rho_dot;
  }

  // Calculate mean predicted measurement
  z_pred.fill(0);
  for (int i = 0; i < Zsig.cols(); i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // Calculate measurement covariance matrix S
  MatrixXd R = MatrixXd(3, 3);
  R.fill(0);
  R(0, 0) = std_radr_ * std_radr_;
  R(1, 1) = std_radphi_ * std_radphi_;
  R(2, 2) = std_radrd_ * std_radrd_;

  S.fill(0);
  for (int i = 0; i < Zsig.cols(); i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // Angle normalization
    while (z_diff(1) >  M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  S = S + R;

  // Calculate cross correlation matrix
  Tc.fill(0);
  for (int i = 0; i < Xsig_pred_.cols(); i++) {
      Tc = Tc + weights_(i) * (Xsig_pred_.col(i) - x_) *
            (Zsig.col(i) - z_pred).transpose();
  }

  // Calculate residual
  VectorXd z_diff = z - z_pred;
  // Angle normalization
  while (z_diff(1) >  M_PI) z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2. * M_PI;

  // Calculate Kalman gain K;
  MatrixXd Sinv = S.inverse();
  MatrixXd K = Tc * Sinv;

  // Update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  // Calculate NIS
  NIS_radar_ = z_diff.transpose() * Sinv * z_diff;
}

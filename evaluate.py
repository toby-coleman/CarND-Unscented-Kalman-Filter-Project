""" Evaluates performance of unscented Kalman filter and plots charts
"""

import pandas as pd
import sys


def read_file(file_name):
    df = pd.read_csv(file_name, delimiter='\t')
    return df

def plot_versus_ground_truth(df):
    import matplotlib
    from matplotlib import pyplot as plt
    plt.ioff()  # Don't use interactive mode

    fig = plt.figure(figsize=(15, 7))
    ax = fig.add_subplot(111)

    # Plot each variable
    ax.plot(df.px_state, df.py_state,
            label='Filter estimates', ls='None', marker='o')
    ax.plot(df.px_ground_truth, df.py_ground_truth,
            label='Ground truth', ls='None', marker='.')
    ax.plot(df.px_measured, df.py_measured,
            label='Measurements', ls='None', marker='x')

    ax.grid()
    ax.set_axis_bgcolor('whitesmoke')
    ax.legend(loc='lower right')
    ax.set_ylabel('y, m')
    ax.set_xlabel('x, m')

    # Save plot
    fig.savefig('figures/filter_output.png')
    plt.close(fig)


def plot_nis_chart(df):
    lidar = df[df.sensor_type == 'lidar']
    radar = df[df.sensor_type == 'radar']

    import matplotlib
    from matplotlib import pyplot as plt
    plt.ioff()  # Don't use interactive mode

    fig = plt.figure(figsize=(15, 7))
    ax = fig.add_subplot(111)

    stats = {}

    if not lidar.empty:
        # Plot lidar NIS
        seconds = (lidar.time_stamp - df.time_stamp.iloc[0]) / 1e6
        line, = ax.plot(seconds, lidar.NIS,
                        label='Lidar', lw=2)
        # Draw line for 5% chi squared critical value
        ax.axhline(5.991, color=line.get_color(), ls='--', lw=3)
        # Compute proportion over 5% level (2 deg freedom)
        stats['lidar'] = len(lidar[lidar.NIS > 5.991]) / len(lidar)
    if not radar.empty:
        seconds = (radar.time_stamp - df.time_stamp.iloc[0]) / 1e6
        line, = ax.plot(seconds, radar.NIS,
                        label='Radar', lw=2)
        # Draw line for 5% chi squared critical value
        ax.axhline(7.815, color=line.get_color(), ls='--', lw=3)
        # Compute proportion over 5% level (3 deg freedom)
        stats['radar'] = len(radar[radar.NIS > 7.815]) / len(radar)

    ax.grid()
    ax.set_axis_bgcolor('whitesmoke')
    ax.legend(loc='upper right')
    ax.set_ylabel('NIS')
    ax.set_ylim([0, 10])
    ax.set_xlabel('Time, seconds')

    # Save plot
    fig.savefig('figures/nis_plot.png')
    plt.close(fig)

    return stats

if __name__ == '__main__':
    # Read in the output file from UKF code
    file_name = sys.argv[1]
    data = read_file(file_name)

    # Plot x and y points
    plot_versus_ground_truth(data)

    # Analyse NIS data
    stats = plot_nis_chart(data)
    print(stats)

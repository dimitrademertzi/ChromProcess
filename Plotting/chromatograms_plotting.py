from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


def plot_peaks_and_boundaries(chroms):
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    z = []
    X, Y = np.meshgrid(chroms[0].time[:12990], np.linspace(0, len(chroms), len(chroms)))
    for count, chrom in enumerate(chroms):
        z.append(chrom.signal[20:12990])

    z = np.array(z)
    ax.plot_wireframe(X, Y, z)
    plt.show()


def heatmap_cluster(chroms, plot_bounds):
    from ChromProcess.Utils.utils.clustering import cluster

    peak_pos = np.array([])
    for chrom in chroms:
        peak_pos = np.hstack((peak_pos, np.array([*chrom.peaks.keys()])))
    peaks = peak_pos[np.argsort(peak_pos)]
    clusters = []
    for c in cluster(peaks, bound=0.02):
        clusters.append(c)

    z = []
    for count, chrom in enumerate(chroms):
        z.append(chrom.signal[plot_bounds[0] : plot_bounds[1]])

    z = np.array(z)
    fig, ax = plt.subplots()
    extent = [
        chroms[0].time[plot_bounds[0]],
        chroms[0].time[plot_bounds[1]],
        0,
        len(chroms),
    ]
    plt.imshow(z, extent=extent, aspect=0.01)
    for c1, clust in enumerate(clusters):
        chrom_plot = []
        rt_plot = []
        for c2, chrom in enumerate(chroms):
            for pk in chrom.peaks.keys():
                if pk in clust:
                    chrom_plot.append(c2)
                    rt_plot.append(pk)
        ax.plot(rt_plot, chrom_plot)

    plt.show()


def peak_area(time, signal, picked_peaks, save_folder=None):
    fig, ax = plt.subplots()
    signal = signal - min(signal)

    ax.plot(time, signal)
    for x in range(0, len(picked_peaks["Peak_indices"])):
        peak_range = range(
            picked_peaks["Peak_start_indices"][x], picked_peaks["Peak_end_indices"][x]
        )
        peak_centre = picked_peaks["Peak_indices"][x]
        ax.fill_between(time[peak_range], signal[peak_range], alpha=0.3)
        ax.plot(
            [time[peak_centre], time[peak_centre]],
            [0, signal[peak_centre]],
            "k--",
            alpha=0.7,
        )
    plt.tight_layout()
    plt.savefig(save_folder)
    plt.close()

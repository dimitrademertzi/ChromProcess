# from pathlib import Path
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from ChromProcess.Loading import conditions_from_csv, mineral_conditions_from_csv, chrom_from_csv
from ChromProcess.Loading.analysis_info.analysis_from_toml import analysis_from_toml
from ChromProcess.Utils.signal_processing.deconvolution import deconvolute_peak
from pathlib import Path
from Plotting.chromatograms_plotting import peak_area
from Plotting.chromatograms_plotting import heatmap_cluster
from ChromProcess.Utils.peak_finding import find_peaks_scipy
from ChromProcess.Utils import indices_from_boundary, peak_indices_to_times
from ChromProcess.Utils import valid_deconvolution
import pandas as pd
from ChromProcess.Processing import add_peaks_to_chromatogram
from ChromProcess.Processing import integrate_chromatogram_peaks
from ChromProcess.Processing import internal_standard_integral_look_ahead
import numpy as np
from ChromProcess import Classes

experiment_number = "MIN001A"
experiment_folder = Path(
    f"{Path.home()}//Macdocs/Master/Internship/Data/{experiment_number}"
)
chromatogram_directory = Path(experiment_folder, f"ChromatogramCSV")
conditions_file = Path(experiment_folder, f"{experiment_number}_conditions.csv")
analysis_file = Path(experiment_folder, f"{experiment_number}_analysis_details.toml")
peak_collection_directory = Path(experiment_folder, f"PeakCollections")

if "FRN" in experiment_number:
    conditions = conditions_from_csv(conditions_file)
else:
    conditions = mineral_conditions_from_csv(conditions_file)

analysis = analysis_from_toml(analysis_file)
if not valid_deconvolution(analysis):
    print("Invalid deconvolution paramaters")
    exit()

# Load in all chromatograms
os.makedirs(peak_collection_directory, exist_ok=True)
chromatogram_files = os.listdir(chromatogram_directory)
chromatogram_files.sort()
chromatogram_files.remove(".DS_Store")
chroms = []
for f in chromatogram_files:
    chroms.append(chrom_from_csv(f"{chromatogram_directory}/{f}"))


# blank = chrom_from_csv(f"{experiment_folder}\\blank_CC.csv")
# subtract baseline
# for chrom in chroms:
#    chrom.signal -= blank.signal


# fig, ax = plt.subplots()
# for c in chroms:
#     ax.plot(c.time[analysis.plot_region[0]:analysis.plot_region[1]],
#            c.signal[analysis.plot_region[0]:analysis.plot_region[1]],
#            label = c.filename)
# plt.show()
# plt.close()

is_start = analysis.internal_standard_region[0]
is_end = analysis.internal_standard_region[1]
for c in chroms:
    c.signal = c.signal - min(
        c.signal[analysis.plot_region[0] : analysis.plot_region[1]]
    )
    internal_standard_integral_look_ahead(c, is_start, is_end)
    c.signal = c.signal / c.internal_standard.height


plot_figures = False
threshold = analysis.peak_pick_threshold
threshold = [threshold for r in analysis.regions]
peak_figure_folder = Path(experiment_folder, "peak_figures")
if type(threshold) == float:
    peak_figure_folder.mkdir(exist_ok=True)
for chrom in chroms:
    for reg, thres in zip(analysis.regions, threshold):
        inds = indices_from_boundary(chrom.time, reg[0], reg[1])
        time = chrom.time[inds]
        signal = chrom.signal[inds]
        picked_peaks = find_peaks_scipy(
            signal,
            threshold=thres,
            min_dist=analysis.peak_distance,
            max_inten=1e100,
            prominence=analysis.prominence,
            wlen=1001,
            look_ahead=analysis.boundary_window,
            smooth_window=11,
        )
        peak_features = peak_indices_to_times(time, picked_peaks)
        peaks = []
        for x in range(0, len(picked_peaks["Peak_indices"])):
            pk_idx = picked_peaks["Peak_indices"][x]
            start_idx = picked_peaks["Peak_start_indices"][x]
            end_idx = picked_peaks["Peak_end_indices"][x]

            retention_time = time[pk_idx]
            start = time[start_idx]
            end = time[end_idx]
            height = signal[pk_idx] - min(
                signal
            )  # subtract the baseline of the region from the peak height
            peaks.append(
                Classes.Peak(retention_time, start, end, indices=[], height=height)
            )
        if plot_figures == True:
            peak_area(
                time,
                signal,
                picked_peaks,
                save_folder=f"{peak_figure_folder}/{reg[0]}_{chrom.filename[:-4]}.png",
            )
        add_peaks_to_chromatogram(peaks, chrom)
    integrate_chromatogram_peaks(chrom, baseline_subtract=True)


### heatmap_cluster(chroms)
for count, reg in enumerate(analysis.deconvolve_regions):
    region_start = analysis.deconvolve_regions[reg]["region_boundaries"][0]
    region_end = analysis.deconvolve_regions[reg]["region_boundaries"][1]
    indices = indices_from_boundary(
        chroms[0].time,
        analysis.deconvolve_regions[reg]["region_boundaries"][0],
        analysis.deconvolve_regions[reg]["region_boundaries"][1],
    )
    peak_folder = f"{experiment_folder}/deconvolved_peaks/{reg, region_start}"
    os.makedirs(peak_folder, exist_ok=True)
    fit_values = np.array(["sample_name", "mse"])
    for n in range(1, analysis.deconvolve_regions[reg]["number_of_peaks"] + 1):
        fit_values = np.hstack((fit_values, [f"amp{n}", f"centre{n}", f"sigma{n}"]))
    fit_values = np.hstack((fit_values, "baseline"))

    # To deconvolve specific samples in a deconvolution region
    chrom_selection = False
    if "selected_chromatograms" in analysis.deconvolve_regions[reg]:
        chrom_selection = True

    for chrom in chroms:
        chrom_filename = chrom.filename.split("_")[1].split(".")[0]
        deconvolve_this = True
        if chrom_selection == True:
            if chrom_filename[:-1] not in analysis.deconvolve_regions[reg]["selected_chromatograms"]:
                deconvolve_this = False
            else:
                exit
    
        if deconvolve_this:
            #filter chromatogram.peaks for any peaks between region start and region end
            #if any peaks is within this margin remove it from chromatogram.peaks
            peaks_to_remove = []
            for peak in chrom.peaks:
                #print(peak)
                if (region_start < peak and  peak < region_end):
                    #print(peak)
                    peaks_to_remove.append(peak)
            for peak in peaks_to_remove:
                chrom.peaks.pop(peak)
    
            popt, pcov, mse, peaks = deconvolute_peak(
                chrom,
                peak_folder,
                indices,
                analysis.deconvolve_regions[reg],
                plotting=True,
            )
            fit_values = np.vstack(
                (
                    fit_values,
                    np.array([chrom.filename.split("_")[1].split(".")[0], mse, *popt]),
                )
            )
            k = [*chrom.peaks.keys()]
            v = [*chrom.peaks.values()]

            for peak in peaks:
                rt = peak.retention_time
                idx = np.where((chrom.time >= peak.start) & (chrom.time <= peak.end))[0]
                peak.indices = idx
                insert = np.searchsorted(k, rt)
                k.insert(insert, rt)
                v.insert(insert, peak)
            chrom.peaks = dict(zip(k, v))

    np.savetxt(
        f"{peak_folder}/gaussian_fit_{region_start}.csv",
        fit_values,
        fmt="%s",
        delimiter=",",
    )


# for chrom in chroms:
#    peaks_indices = peak_indices_from_file(chrom,f"{peak_collection_directory}\\{chrom.filename}")
#    peak_starts, peak_ends = peak_boundaries_from_file(chrom,f"{peak_collection_directory}\\{chrom.filename}")
#    picked_peaks = {'Peak_indices':peaks_indices, 'Peak_start_indices':peak_starts, 'Peak_end_indices':peak_ends}
#
#    peak_features = peak_indices_to_times(chrom.time,picked_peaks)
#    add_peaks_to_chromatogram(peak_features, chrom)
#    integrate_chromatogram_peaks(chrom)

colors = []
color_palette_list = [
    "Reds",
    "Purples",
    "RdPu",
    "Wistia",
    "Blues",
    "YlGn",
    "gray_r",
]  # "copper_r", "GnBu", "BuPu", , "hls", "Set2", "RdPu", ]
for color in color_palette_list:
    colors += sns.color_palette(f"{color}", 5).as_hex()
colors2 = colors[::-1]

# sns.set_style("dark")
fig, ax = plt.subplots()
ax.set_prop_cycle(color=[c for c in colors2])
for c in chroms:
    ax.plot(
        c.time[analysis.plot_region[0] : analysis.plot_region[1]],
        c.signal[analysis.plot_region[0] : analysis.plot_region[1]],
        label=c.filename.split("_")[1].split(".")[0],
    )
handles, labels = ax.get_legend_handles_labels()
ax.legend(
    handles, labels, ncol=2, fontsize=10, bbox_to_anchor=(1.1, 1.1), loc="upper right"
)
# ax.set_xlim(12.0, 13.0)
# ax.set_ylim(0, 0.5) #A series (-0.05), B series (-0.0025), C series (0)
plt.show()

# heatmap_cluster(chroms,analysis.plot_region, peak_agglomeration_boundary=0.02)
#for c, v in zip(chroms, conditions.series_values):
#    c.write_peak_collection(
#        filename=f"{peak_collection_directory}/{c.filename}",
#        header_text=f"{conditions.series_unit},{v}\n",
#    )

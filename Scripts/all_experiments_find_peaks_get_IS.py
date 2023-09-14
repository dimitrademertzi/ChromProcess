# from pathlib import Path
import os
import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from ChromProcess.Loading import (
    conditions_from_csv,
    mineral_conditions_from_csv,
    chrom_from_csv,

    peak_collection_from_csv,
    mineral_peak_collection_from_csv,
)
from ChromProcess.Loading.analysis_info.analysis_from_toml import analysis_from_toml
from pathlib import Path
from Plotting.chromatograms_plotting import peak_area
from ChromProcess.Utils.peak_finding import find_peaks_scipy
from ChromProcess.Utils import indices_from_boundary, peak_indices_to_times
import pandas as pd
from ChromProcess.Processing import add_peaks_to_chromatogram
from ChromProcess.Processing import integrate_chromatogram_peaks
from ChromProcess.Processing import internal_standard_integral_look_ahead
import numpy as np
from ChromProcess import Classes
from functools import reduce
import seaborn as sns

peak_tables = []
sample_names = []
chroms: list[Classes.Chromatogram] = []
experiment_number = ["MIN003", "MIN004", "MIN005", "MIN008", "MIN009", "MIN010", "MIN011", 'PIN001']


peak_borders_list = [9.75, 10.05, 
                     11.15, 11.35, 
                     11.5, 11.65, 
                     11.75, 12.12,
                     13.55, 13.67,
                     14.02, 15.05,
                     15.52, 15.85,
                     16.5, 16.84,
                     17.07, 17.42,
                     17.55, 18.00,
                     19.7, 19.77]

colors_single = ["Reds", "Greens", "magma", "Purples", "Blues", "Wistia", "gray_r", "Oranges"]
colors_list = []
colors_final  = []

for actual_number, number in enumerate(experiment_number):
    experiment_folder = Path(
        f"{Path.home()}//Macdocs/Master/Internship/Data/{number}"
    )
    chromatogram_directory = Path(experiment_folder, f"ChromatogramCSV")
    conditions_file = Path(experiment_folder, f"{number}_conditions.csv")
    analysis_file = Path(experiment_folder, f"{number}_analysis_details.toml")
    data_report_directory = Path(f"{experiment_folder}/DataReports")
    peak_collection_directory = Path(f"{experiment_folder}/PeakCollections2")
    
    if "FRN" in experiment_number:
        conditions = conditions_from_csv(conditions_file)
    else:
        conditions = mineral_conditions_from_csv(conditions_file)

    analysis = analysis_from_toml(analysis_file)
    chromatogram_files = os.listdir(chromatogram_directory)
    chromatogram_files.sort()

    if ".DS_Store" in chromatogram_files:
        chromatogram_files.remove(".DS_Store")

    for n, f in enumerate(chromatogram_files):
        chroms.append(chrom_from_csv(f"{chromatogram_directory}/{f}"))
        colors_list.append(colors_single[actual_number])
    
    is_start = analysis.internal_standard_region[0]
    is_end = analysis.internal_standard_region[1]


            
    for chrom in chroms:
        chrom.experiment_name = chrom.filename.split(".")[0:2]
        chrom.signal = chrom.signal - min(
        chrom.signal[analysis.plot_region[0] : analysis.plot_region[1]])
        

color_dict = {'MIN003':'#440154', 
              'MIN004':'#46327e', 
              'MIN005':'#365c8d', 
              'MIN006':'#277f8e', 
              'MIN008':'#1fa187', 
              'MIN009':'#4ac16d',
              'MIN010':'#a0da39',
              'MIN011':'#fde725',
              'PIN001':'#440154'
              }

fig, ax = plt.subplots()
for color in colors_list:
    colors_final +=  sns.color_palette(color, 1).as_hex()
ax.set_prop_cycle(color=[c for c in colors_final])

exp_list = ['MIN009', 'MIN010', 'MIN011', 'PIN001']
for c in chroms:
    if '(' in c.experiment_name:
        exit
    else:
        if c.experiment_name[0] in exp_list:
            ax.plot(
                c.time[analysis.plot_region[0] : analysis.plot_region[1]],
                c.signal[analysis.plot_region[0] : analysis.plot_region[1]],
                label=c.experiment_name[0], color=color_dict[c.experiment_name[0]]
            )
            handles, labels = ax.get_legend_handles_labels()
            unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
            ax.legend(
                *zip(*unique), ncol=2, fontsize=10, bbox_to_anchor=(1.1, 1.1), loc="upper right"
            )
plt.vlines(x = peak_borders_list, ymin=-10, ymax=400, color='k')
plt.show()

 
#IS_heights = []
#IS_integrals = []
#experiment_number_list = []
#IS_df = pd.DataFrame(columns=('experiment_number', 'IS_height', 'IS_integral'))
#
#experiment_number = ["MIN003", "MIN004", "MIN005", "MIN006", "MIN008", "MIN009"]
#peak_files_all = []
#
#for actual_number, number in enumerate(experiment_number):
#    experiment_folder = Path(
#        f"{Path.home()}//Macdocs/Master/Internship/Data/{number}"
#    )
#    chromatogram_directory = Path(experiment_folder, f"ChromatogramCSV")
#    conditions_file = Path(experiment_folder, f"{number}_conditions.csv")
#    analysis_file = Path(experiment_folder, f"{number}_analysis_details.toml")
#    data_report_directory = Path(f"{experiment_folder}/DataReports")
#    peak_collection_directory = Path(f"{experiment_folder}/PeakCollections")
#    
#    os.makedirs(data_report_directory, exist_ok=True)
#    peak_files = os.listdir(peak_collection_directory)
#
#    if "FRN" in experiment_number:
#        conditions = conditions_from_csv(conditions_file)
#    else:
#        conditions = mineral_conditions_from_csv(conditions_file)
#
#    analysis = analysis_from_toml(analysis_file)
#    peak_tables = []
#    peak_files = os.listdir(peak_collection_directory)
#
#
#    if ".DS_Store" in peak_files:
#        peak_files.pop(peak_files.index(".DS_Store"))
#    peak_files_all.append(peak_files) 
#
#    for file in peak_files:
#        if file.endswith(".csv") or file.endswith(".CSV"):
#            if "FRN" in experiment_number:
#                peak_tables.append(
#                    peak_collection_from_csv(
#                        f"{peak_collection_directory}/{file}", round_digits=7
#                    )
#                )
#            else:
#                peak_tables.append(
#                    mineral_peak_collection_from_csv(
#                    f"{peak_collection_directory}/{file}", round_digits=7
#                    )
#                )
#
#    series = Classes.PeakCollectionSeries(
#        peak_tables, name=f"{experiment_number[actual_number]}", conditions=conditions.conditions
#    )
#
#    
#    IS_heights.append([i.internal_standard.height for i in series.peak_collections])
#    IS_integrals.append([i.internal_standard.integral for i in series.peak_collections])
#
#
#flattened_peak_files = reduce(lambda z, y: z + y, peak_files_all)
#experiment_number_list.append([filename[:6] for filename in flattened_peak_files])
#
#flattened_exp_num, flattened_IS_height, flattened_IS_integral  = reduce(lambda z, y: z + y, experiment_number_list), reduce(lambda z, y: z + y, IS_heights), reduce(lambda z, y: z + y, IS_integrals)
#
#IS_df[['experiment_number', 'IS_height', 'IS_integral']] = list(zip(flattened_exp_num, flattened_IS_height, flattened_IS_integral))
#
#fig, ax = plt.subplots(2, 1, sharex=True)
#
#for count, column_name in enumerate(IS_df.columns):
#    if column_name != 'experiment_number':
#        sns.set_palette('viridis')
#        sns.barplot(data=IS_df, x=IS_df.experiment_number, y=IS_df[column_name], ax=ax[count-1])
#        ax[0].set(xlabel=None)
#        ax[0].tick_params(axis='x', which='both', bottom=False)
#
#
#plt.show()

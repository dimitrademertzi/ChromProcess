import os
from ChromProcess import Classes
from ChromProcess.Loading import (
    peak_collection_from_csv,
    mineral_peak_collection_from_csv,
)
from ChromProcess.Loading import analysis_from_toml
from ChromProcess.Loading import conditions_from_csv, mineral_conditions_from_csv
from pathlib import Path

experiment_number = "MIN008"
experiment_folder = Path(
    f"{Path.home()}//Macdocs/Master/Internship/Data/{experiment_number}"
)
peak_collection_directory = Path(f"{experiment_folder}/PeakCollections")
conditions_file = Path(f"{experiment_folder}/{experiment_number}_conditions.csv")
analysis_file = Path(f"{experiment_folder}/{experiment_number}_analysis_details.toml")
data_report_directory = Path(f"{experiment_folder}/DataReports")
os.makedirs(data_report_directory, exist_ok=True)

if "FRN" in experiment_number:
    conditions = conditions_from_csv(conditions_file)
else:
    conditions = mineral_conditions_from_csv(conditions_file)
analysis = analysis_from_toml(analysis_file)
peak_tables = []
peak_files = os.listdir(peak_collection_directory)

if ".DS_Store" in peak_files:
    peak_files.pop(peak_files.index(".DS_Store"))

for file in peak_files:
    if file.endswith(".csv") or file.endswith(".CSV"):
        if "FRN" in experiment_number:
            peak_tables.append(
                peak_collection_from_csv(
                    f"{peak_collection_directory}/{file}", round_digits=7
                )
            )
        else:
            peak_tables.append(
                mineral_peak_collection_from_csv(
                f"{peak_collection_directory}/{file}", round_digits=7
                )
            )

# Create series of peak collections
series = Classes.PeakCollectionSeries(
    peak_tables, name=f"{experiment_number}", conditions=conditions.conditions
)

IS_pos = 9.33

# series.reference_integrals_to_IS()
# 5% of internal standard integral if integrals are normalised to IS

peak_removal_limit = 0.001

series.remove_peaks_below_threshold(peak_removal_limit, metric="height")
peak_agglomeration_boundary = 0.015  # distance cutoff
# cluster_threshold = 0.008

series.get_peak_clusters(bound=peak_agglomeration_boundary)
# to_remove = []
# for c1, clust in enumerate(series.clusters):
#    max_height = 0
#    for pc in series.peak_collections:
#        for pk in pc.peaks:
#            if pk.retention_time in clust and pk.height>max_height:
#                max_height = pk.height
#        #for pk in pc.peaks:
#        #    if pk.height < (max_height/50):
#
#    if max_height < cluster_threshold:
#        to_remove.append(c1)
#
# [series.clusters.pop(c) for c in sorted(to_remove,reverse=True)]
#        to_remove = []
#        for k in integral_dict:
#            if max(integral_dict[k]) < cluster_removal_limit:
#                to_remove = to_remove + [k]
#        [integral_dict.pop(key) for key in to_remove]

if "FRN" in experiment_number:
    series.write_data_reports(
        f"{data_report_directory}/{series.name}", analysis
    )  # create arrays for output
else:
    series.mineral_write_data_reports(
        f"{data_report_directory}/{series.name}", analysis
    )  # create arrays for output

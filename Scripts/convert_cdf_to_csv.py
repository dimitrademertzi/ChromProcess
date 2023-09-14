# from pathlib import Path
import os
import matplotlib.pyplot as plt
from ChromProcess.Loading import conditions_from_csv, chrom_from_csv, chrom_from_cdf
from ChromProcess.Loading.analysis_info.analysis_from_toml import analysis_from_toml
from ChromProcess.Processing import ic_background_subtraction
from ChromProcess.Writers import chromatogram_to_csv
from pathlib import Path
import pandas as pd
import numpy as np



experiment_number = "MIN012"
experiment_folder = Path(
    f"{Path.home()}//Macdocs/Master/Internship/Data/{experiment_number}"
)
chromatogram_directory = Path(experiment_folder, f"CDF_files")
csv_directory = Path(experiment_folder, "ChromatogramCSV_GCMS")
os.makedirs(csv_directory, exist_ok=True)

analysis_file = Path(experiment_folder, f"{experiment_number}_analysis_details.toml")
analysis = analysis_from_toml(analysis_file)

chromatogram_files = os.listdir(chromatogram_directory)
chromatogram_files.sort()
#chromatogram_files.remove(".DS_Store")
chroms = []

for f in chromatogram_files:
    chroms.append(chrom_from_cdf(f"{chromatogram_directory}/{f}", load_ms=True))

for n,c in enumerate(chroms):
    chroms[n].signal = ic_background_subtraction(chromatogram=chroms[n], threshold=50)
    chromatogram_to_csv(chroms[n], f"{csv_directory}/{chroms[n].filename[:-4]}.csv")
    print(chroms[n].filename)


fig, ax = plt.subplots()
for n, c in enumerate(chroms):
    if True:
        print(c.filename)
        ax.plot(
                    c.time,
                    c.signal,
                    label=c.filename[:-4]
                )
        
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles = handles, labels=labels)
#plt.show()
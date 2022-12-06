
# from pathlib import Path
import os
import matplotlib.pyplot as plt
from ChromProcess.Loading import conditions_from_csv, chrom_from_csv, chrom_from_cdf
from ChromProcess.Loading.analysis_info.analysis_from_toml import analysis_from_toml
from ChromProcess.Utils.signal_processing.deconvolution import deconvolute_peak
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

experiment_number = "FRN154"
experiment_folder = Path(f"{Path.home()}/PhD/Data/{experiment_number}")
chromatogram_directory = Path(experiment_folder, f"cdf")
csv_directory = Path(experiment_folder, "chromatogramCSV")
conditions_file = Path(experiment_folder, f"{experiment_number}_conditions.csv")
analysis_file = Path(experiment_folder, f"{experiment_number}_analysis_details.toml")
peak_collection_directory = Path(experiment_folder, f"PeakCollections")

conditions = conditions_from_csv(conditions_file)
analysis = analysis_from_toml(analysis_file)

os.makedirs(csv_directory, exist_ok=True)
chromatogram_files = os.listdir(chromatogram_directory)
chromatogram_files.sort()
chroms = []
for f in chromatogram_files:
    chroms.append(chrom_from_cdf(f"{chromatogram_directory}/{f}"))

fig, ax = plt.subplots()
for c in chroms:
    ax.plot(
        c.time[analysis.plot_region[0] : analysis.plot_region[1]],
        c.signal[analysis.plot_region[0] : analysis.plot_region[1]],
        label=c.filename,
    )
plt.show()

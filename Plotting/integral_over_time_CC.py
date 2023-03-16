import matplotlib.pyplot as plt
import numpy as np
from ChromProcess.Loading.data_report import data_report_from_csv

data_folder = "../../../Data/"
experiment_list = ["CC01", "CC02", "CC03", "CC04", "CC05", "CC06", "CC07", "CC08"]
for exp in experiment_list:
    exp_folder = f"{data_folder}{exp}"
    data_report = data_report_from_csv(
        f"{exp_folder}/DataReports/{exp}_HPLC_integral_report.csv"
    )
    peptide = data_report.experiment_code
    data = data_report.data
    fig, ax = plt.subplots(dpi=300)
    times = data_report.series_values
    for peak in data:
        normalized_integrals = data[peak] / max(data[peak])
        ax.plot(times, normalized_integrals, "o--", label=peak, alpha=0.8)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Normalized integrals")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"{data_folder}/CC01/{peptide}_integral_plot.png")
    plt.close()

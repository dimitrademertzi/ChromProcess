import os
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
from ChromProcess import Classes
from os.path import isfile, join
from ChromProcess.Loading import chrom_from_csv
from ChromProcess.Processing import internal_standard_integral_look_ahead
from ChromProcess.Loading.analysis_info.analysis_from_toml import analysis_from_toml

color_dict = {
    "MIN003": "#440154",
    "MIN004": "#46327e",
    "MIN005": "#365c8d",
    "MIN006": "#277f8e",
    "MIN008": "#1fa187",
    "MIN009": "#4ac16d",
    "MIN010": "#a0da39",
    "MIN011": "#fde725",
    "PIN001": "#440154",
}

experiments = ["003", "004", "005", "006", "008", "009", "010", "011"]
chroms: list[Classes.Chromatogram] = []

for exp in experiments:
    exp_dir = Path(os.environ["INTERNSHIP"]) / "Data" / f"MIN{exp}"
    chroms_dir = (
        Path(os.environ["INTERNSHIP"]) / "Data" / f"MIN{exp}" / "ChromatogramCSV"
    )

    analysis_file = Path(exp_dir, f"MIN{exp}_analysis_details.toml")
    analysis = analysis_from_toml(analysis_file)

    chrom_files = os.listdir(chroms_dir)
    chrom_files.sort()

    if ".DS_Store" in chrom_files:
        chrom_files.remove(".DS_Store")

    for f in chrom_files:
        chroms.append(chrom_from_csv(f"{chroms_dir}/{f}"))

    for c in chroms:
        c.mineral = c.filename.split("_")[1]
        c.experiment_name = c.filename.split("_")[0]

        c.singal = c.signal - min(
            c.signal[analysis.plot_region[0] : analysis.plot_region[1]]
        )
        # c.signal = c.signal / c.internal_standard.height

exp_list = []
sns.set_style("dark")
fig, ax = plt.subplots()


peak_rt = 20
boundary = 0.05

for c in chroms:
    if c.mineral not in exp_list:
        ax.plot(
            c.time[analysis.plot_region[0] : analysis.plot_region[1]],
            c.signal[analysis.plot_region[0] : analysis.plot_region[1]],
            color=color_dict[c.experiment_name],
            linestyle="solid" if "(" not in c.mineral else "None",
        )
        ax.set_xlim(peak_rt - boundary, peak_rt + boundary)
plt.show()

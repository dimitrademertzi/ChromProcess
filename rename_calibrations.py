from pathlib import Path
import pandas as pd
import csv
import shutil

number_to_sugar = {
    "1": "lyxose",
    "2": "xylose",
    "3": "glycerol",
    "4": "arabitol",
    "5": "dulcitol",
    "6": "sorbitol",
    "7": "xylitol",
    "8": "ribitol",
    "9": "erythritol",
    "10": "maltitol",
}

initial_concentrations = {
    "lyxose": 2.0155,
    "xylose": 1.9875,
    "glycerol": 2.1955,
    "arabitol": 2.0225,
    "dulcitol": 2.003,
    "sorbitol": 1.9915,
    "xylitol": 2.035,
    "ribitol": 1.9865,
    "erythritol": 2.0365,
    "maltitol": 2.0025,
}

final_concentrations = {
    "lyxose": {},
    "xylose": {},
    "glycerol": {},
    "arabitol": {},
    "dulcitol": {},
    "sorbitol": {},
    "xylitol": {},
    "ribitol": {},
    "erythritol": {},
    "maltitol": {},
}

root = Path(
    f"{Path.home()}/Macdocs/Master/Internship/Data/Calibrations/ChromatogramCSV"
)
new_folders = Path(f"{Path.home()}/Macdocs/Master/Internship/Data/Calibrations")
assert root.exists()
files = root.glob("*.CSV")
for f in files:
    name = f.name.split(".")[0]
    if not ("blank" in name or "kovats" in name):
        if "CAL" in name:
            name = name.split("_")[1]
        conv_name = f"{number_to_sugar[name[:-2]]}_{name[-2]}_{name[-1]}.csv"
        subfolder = Path(f"{new_folders}/{number_to_sugar[name[:-2]]}/chromatograms/")
        subfolder.mkdir(exist_ok=True, parents=True)
        shutil.copy(
            f, f"{new_folders}/{number_to_sugar[name[:-2]]}/chromatograms/{conv_name}"
        )

        concentration = initial_concentrations[number_to_sugar[name[:-2]]] / (
            (2 ** int(name[-1])) * 1000
        )
        concentration = str(concentration)[:10]
        final_concentrations[number_to_sugar[name[:-2]]][conv_name] = concentration

for key in final_concentrations.keys():
    name = f.name.split(".")[0]
    print(key)

    with open(Path(f"{new_folders}/{key}/{key}_concentrations.csv"), "w") as conc_file:
        writer = csv.writer(conc_file)
        header = ("filename", r"concentration/ M")
        writer.writerow(header)
        writer.writerows(final_concentrations[key].items())

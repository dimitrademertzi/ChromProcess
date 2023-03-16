import seaborn as sns

colors = []
color_palette_list = [
    "Greys",
    "YlGn",
    "YlGnBu",
    "gist_heat_r",
    "magma_r",
    "mako_r",
    "ocean",
    "rocket_r",
    "PuRd",
    "BuPu",
    "hot_r",
    "autumn_r",
    "YlOrRd",
    "GnBu",
]
for color in color_palette_list:
    colors += sns.color_palette(f"{color}", 5).as_hex()

# for c in colors:
#    print(c)

current_color_codes = sns.color_palette("hot_r", 8).as_hex()
print(current_color_codes)

import numpy as np
import re
import matplotlib.pyplot as plt
from data_prep import data
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

data_dict = data.get_prepared_data(x_time_path="/home/william/Desktop/Uni/ppn/nuppn/frames/ppn/i_process/x-time.dat")
# data_dict["SM 144"]["abundance"] = ["1"]
# print(data_dict["SM 144"])

print(data_dict["PR 141"])
# input()

def split(elem):
    reg = r"^(?P<elem>[a-zA-Z]+)\s*(?P<mass>\d+)"
    reg = re.compile(reg)
    e = re.fullmatch(reg, elem)
    if e is None:
        return None
    return [e['elem'], e['mass']]




def plot_thing(data_dict, neut_range, prot_range):
    stable_prot, stable_neut, stable_name = [], [], []
    isotope_prot, isotope_neut, isotope_name, isotope_abundance = [], [], [], []
    stable_abund = []
    for x, i in enumerate(data_dict):
        if x > 6:
            if i != "NEUT" and i != "PROT":
                if data_dict[i]["sol_abundance"] is not None:
                    stable_prot.append(data_dict[i]["proton"])
                    stable_neut.append(data_dict[i]["neutron"])
                    try:
                        stable_name.append(split(i)[0])
                    except:
                        print(i)
                    stable_abund.append(data_dict[i]["abundance"])
                if data_dict[i]["abundance"] is not None and len(data_dict[i]["abundance"]) != 0:
                    isotope_prot.append(data_dict[i]["proton"])
                    isotope_neut.append(data_dict[i]["neutron"])
                    isotope_name.append(split(i)[0])
                    isotope_abundance.append(data_dict[i]["abundance"][1]+1e-99)


    fig, ax = plt.subplots(figsize=(10, 10))

    def nuclide_plotting(protons, neutrons, names, stable, abundance):
        square = []
        invalid = []


        for q, i in enumerate(protons):
            if i is None:
                invalid.append(q)
        for i in sorted(invalid, reverse=True):
            del protons[i]
            del neutrons[i]
            del names[i]
            del abundance[i]
        for n, z in zip(neutrons, protons):
            if n is not None and z is not None:
                rect = Rectangle((n, z), 1, 1, fill=None)
                square.append(rect)


        abund = []
        for i in range(len(square)):
            rx, ry = square[i].get_xy()
            cx = rx + square[i].get_width() / 2.0
            cy = ry + square[i].get_width() / 2.0
            # if stable is not True:
            plot_name = "$^{" + str(neutrons[i]+protons[i]) + "}_{" + str(protons[i]) +"}{" + names[i] + "}$"
            abund.append(abundance[i])
            # ax.annotate(plot_name, (cx, cy), color="blue", weight="bold", ha="center", va="center", fontsize=2)
        if stable:
            patches = PatchCollection(square, edgecolors="black", match_original=True, lw=3)
        else:
            patches = PatchCollection(square, edgecolors="blue", match_original=True, lw=1, cmap=plt.cm.hot_r)
            patches.set_clim((-30, 0))
            patches.set_array(np.log10((abund)))
            cbar = plt.colorbar(patches, orientation="vertical")
            cbar.set_label("log X$_\mathrm{i}$")
        return patches, cx, cy

    stable_plot_patches, cx, cy = nuclide_plotting(stable_prot, stable_neut, stable_name, True, stable_abund)
    isotope_plot_patches, cx, cy = nuclide_plotting(isotope_prot, isotope_neut, isotope_name, False, isotope_abundance)
    ax.set_ylim(prot_range)
    ax.set_xlim(neut_range)
    ax.add_collection(isotope_plot_patches)
    ax.add_collection(stable_plot_patches)



    return


plot_thing(data_dict, (0, 60), (0, 60))
plt.show()
print("done")

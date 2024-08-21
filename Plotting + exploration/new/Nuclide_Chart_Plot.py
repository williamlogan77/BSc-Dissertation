import numpy as np
import matplotlib.pyplot as plt
import rjs_ppn as ppn
import os
import re
from tqdm import tqdm

temperature, density = "3e8", "1e3"
ele = "Sm"
my_data, my_headers = ppn.read_xtime(f"../ppn/nuppn/frames/ppn/i_process_runs/t{temperature}_p{density}/x-time.dat")

# ~~~~~~~~~~ ABUNDANCE PLOTTING ~~~~~~~~~~

# plt.figure(1)
# isotopes = ["n"]
# ppn.plot_abundance_time(isotopes, my_data, my_headers)
# plt.ylim(1e-14, 1)
# plt.show()


# '''
# ~~~~~~~~~~ NUCLIDE CHART PLOTTING ~~~~~~~~~~~~~
# dir = f"../ppn/nuppn/frames/ppn/i_process_runs/t{temperature}_p{density}/"
# save_dir = f"./plots/i_process/t{temperature}_p{density}/{ele}/"
# for filename in tqdm(os.listdir(dir)):
#     plt.figure(1)
#     if filename.startswith(r"flux"):
#         index = int(re.findall(r"flux_(\d{5})", filename)[0])
#         savefilename = filename.split(".")[0] + ".jpg"
#         # if index <= 40:
#         try:
#             if not os.path.exists(save_dir + savefilename):
#                 # print(save_dir + savefilename, os.path.exists(save_dir + savefilename))
#                 file = dir + "/" + filename
#                 flux = np.genfromtxt(file, skip_header=1)
#                 plt.figure(1)
#                 ppn.plot_nuclide_chart(my_data, index, my_headers, flux=flux, xlim=(84, 96), ylim=(60, 64))
#                 plt.title(f"{filename.split('.')[0]}")
#                 plt.show()
#                 # plt.savefig(f"{save_dir}{filename.split('.')[0]}.jpg")
#                 print('\n' + f"done {filename}")
#                 plt.clf()
#         except TypeError as e:
#             print(e)
#             print(f"this {filename} is bad")

# '''
# ordered_headers = ppn.reorder_isotopes(my_data, my_headers, 1)[1]
# ordered_data = ppn.reorder_isotopes(my_data, my_headers, 1)[0]
# ele_to_search = "Pm-154"

# def plot_relative_abundance(data, headers, element_wanted):
#     for x, i in enumerate(headers):
#         if i == element_wanted:
#             ele_index, name = x, i
#     return(ele_index, name)
#
# print(plot_relative_abundance(ordered_data, ordered_headers, ele_to_search))

# ppn.plot_relative_abundance


fig = plt.subplots(figsize=(15, 15))

input_data = [my_data, my_data, my_data]
input_headers = [my_headers, my_headers, my_headers]
cycles = ["200", "250", "300"]
ppn.plot_elements(input_data, input_headers, cycles, plot="[A/B]", reference="Fe", decay=[False], ylim=(-5, 8))
plt.show()

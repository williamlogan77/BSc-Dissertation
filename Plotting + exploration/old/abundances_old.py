from data_prep import data
import numpy as np
import matplotlib.pyplot as plt
import re
import time
import seaborn as sns

start_time = time.time()

data_dict = data.get_prepared_data('../Data_files/x-time.dat', flux_path='../Data_files/flux_00010.DAT')

def sum_abundances(element, cycle_num=0):
    search_element = element.split()[0]
    same_element = []
    abundance_total = []
    for x in data_dict:
        found_element = x.split()[0]
        if search_element == found_element:
            same_element.append(x)
            if len(data_dict[x]["abundance"]) != 0:
                abundance_total.append(data_dict[x]["abundance"][cycle_num])
    return same_element, np.sum(abundance_total)


def plot_elements(data_dict, e_to_plot):
    fig = plt.figure(figsize=(18, 10))
    x, y = [], []
    for element in e_to_plot:
        element_names, elemental_abundance = sum_abundances(element)
        elemental_abundance = np.log10(elemental_abundance)
        x.append(data_dict[element]["proton"])
        y.append(elemental_abundance)
        plt.scatter(data_dict[element]["proton"], elemental_abundance, label=split(element_names[0])[0])
    plt.plot(x, y)

    plt.ylim([-10, 0])
    plt.tight_layout()
    plt.legend(loc="upper right")
    plt.show()

def plot_abundances(data_dict, e_to_plot):
    fig = plt.figure(figsize=(10, 10))
    for element in e_to_plot:
        if len(data_dict[element]["abundance"]) != 0:
            print(element)
            print(data_dict[element])
            plt.plot(np.log10(data_dict["time"]), np.log10(data_dict[element]["abundance"]), label=element)

    # plt.ylim([-10, 0])
    plt.tight_layout()
    plt.legend(loc="upper right")
    plt.show()



def split(elem):
    reg = r"^(?P<elem>[a-zA-Z]+)\s*(?P<mass>\d+)"
    reg = re.compile(reg)
    e = re.fullmatch(reg, elem)
    if e is None:
        return None
    return [e['elem'], e['mass']]


def get_elements_by_proton(data_dict, p_range):
    """
    :param data_dict: data_dict
    :param range: Tuple (start, end)
    :return:
    """

    def find_by_proton_number(row, proton_number):
        if type(data_dict[row]) != dict:
            return False
        if data_dict[row]["proton"] is None:
            return False
        return data_dict[row]["proton"] == proton_number
    elements_in_range = []
    for i in range(*p_range):
        found_elements = list(filter(lambda x: find_by_proton_number(x, i), data_dict))
        if len(found_elements) >= 1:
            elements_in_range.append(found_elements[0])
    return elements_in_range


elements = ["NEUT"]
plot_abundances(data_dict, elements)


# protons = (24, 36)
# # elements_to_plot = get_elements_by_proton(data_dict, protons)
# elements_to_plot = ["C 12", "O 16", "F 18"]
# plot_abundances(data_dict, elements_to_plot)
#
# #elements_to_plot += get_elements_by_proton(data_dict, (20, 30))
# print(elements_to_plot)
# #elements_to_plot.remove("BE 7")
# plot_elements(data_dict, elements_to_plot)
#
# print(data_dict["HE 3"]["flux"])
# print(data_dict["HE 4"]["flux"])

# def heatmap(data_dict, graphlimits):
#     flux_plot = []
#     proton = []
#     neutron = []
#     for x, i in enumerate(data_dict):
#         if x > 7:
#             if graphlimits[1] >= x >= graphlimits[0]:
#                 flux_plot.append(data_dict[i]["flux"])
#                 proton.append(data_dict[i]["proton"])
#                 neutron.append(data_dict[i]["neutron"])
#     # fig = plt.figure()
#     # plt.imshow(flux_plot)
#     # plt.show()
#     return(flux_plot, proton, neutron)
#
#
# print(heatmap(data_dict, (0, 10)))





print(time.time() - start_time)
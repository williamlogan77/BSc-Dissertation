import numpy as np


def find_by_mass(data_dict, mass, prot):
    for elem in data_dict.keys():
        if type(data_dict[elem]) == np.ndarray:
            continue
        if (data_dict[elem]["mass"] == mass) and (data_dict[elem]["proton"] == prot):
            return elem
    return False


def unpack_flux(lines, data_dict):
    # Z is proton
    # A is total mass
    for x, i in enumerate(lines):
        if x > 1:
            Z, A = i[9:16].split()
            Z, A = int(Z), int(A)
            flux = np.float64(i[62:73])
            elem = find_by_mass(data_dict, A, Z)
            if elem is not False:
                data_dict[elem]["flux"] = flux
    return data_dict
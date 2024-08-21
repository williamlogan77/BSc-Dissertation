import re
import numpy as np

not_elements = ['cycle', 'time', 't9', 'rho', '1-sum(yps)', 'ye', 'NEUT', 'PROT']


def fill_blanks(data):
    def missing_proton(row):
        if type(data[row]) == np.ndarray:
            return False
        return data[row]["proton"] is None

    unknown_proton = list(filter(missing_proton, data))
    all_elements = data.keys() - not_elements  # remove the not elements from the list of actual elements
    for elem in unknown_proton:
        reg = r"^(?P<elem>[a-zA-Z]+)\s*(?P<mass>\d+)"
        reg = re.compile(reg)
        e = re.fullmatch(reg, elem)
        if e is None:
            print(f"{elem} can't be matched with regex")
            continue
        base_elem = e['elem']
        total_mass = int(e['mass'])
        r = re.compile(f"^{base_elem} \d+")
        new_list = list(filter(r.match, all_elements))
        if elem in new_list:
            new_list.remove(elem)
        for element in new_list:
            if data[element]["proton"] is not None:
                data[elem]["proton"] = int(data[element]["proton"])
                data[elem]["neutron"] = total_mass - int(data[element]["proton"])
                data[elem]["mass"] = total_mass
                break

    return data

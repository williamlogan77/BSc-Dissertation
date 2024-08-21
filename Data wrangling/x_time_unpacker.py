import numpy as np
import re

not_elements = ['cycle', 'time', 't9', 'rho', '1-sum(yps)', 'ye']

def get_data_from_lines(lines):
    """ Header part """
    header = lines[0]
    header_elems = header.split("|")
    actual_names = []
    for i, name in enumerate(header_elems[1:]):
        t_name = ' '.join(name.split())
        if i > 5:
            t_name = re.sub(r'^\d+-', '', t_name)
        actual_names.append(t_name)

    """ Data part """
    data = []
    for line in lines[1:]:
        data.append(np.float64(line.split()))
    data = np.array(data)

    """ Dictionary part """
    abundances = {}
    for name, measurements in zip(actual_names, np.transpose(data)):
        try:  # "cycle" doesn't have a number :)
            total = int(name.split()[-1])
        except ValueError:
            total = None
        if name in not_elements:
            abundances[name] = measurements
        else:
            abundances[name] = {
                'proton': None,
                'neutron': None,
                'mass': total,
                'sol_abundance': None,
                'abundance': np.float64(measurements)
            }
    return abundances

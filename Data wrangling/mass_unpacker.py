import numpy as np


def unpack_mass(lines, x_time):
    # exception for the "PROT" line:
    stupid_line = lines[0].split()
    x_time["NEUT"] = {
        'proton': 0,
        'neutron': 1,
        'mass': 1,
        'sol_abundance': 0,
        'abundance': x_time["NEUT"]
    }
    x_time["PROT"] = {
        'proton': 1,
        'neutron': 0,
        'mass': 1,
        'sol_abundance': np.float64(stupid_line[2]),
        'abundance': x_time["PROT"]
    }



    # Normal format:
    for line in lines[1:]:
        p_num = line[:3]
        name = line[3:6]
        total = line[6:9]
        n_num = int(total) - int(p_num)
        actual_name = name.strip().upper() + ' ' + total.strip()
        if actual_name in x_time.keys():
            data = x_time[actual_name]['abundance']
        else:
            data = []
        x_time[actual_name] = {
            'proton': int(p_num),
            'neutron': int(n_num),
            'mass': int(total),
            'sol_abundance': np.float64(line.split()[-1]),
            'abundance': np.float64(data)
        }


    return x_time

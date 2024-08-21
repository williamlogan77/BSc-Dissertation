import numpy as np
import matplotlib as plt
import math
import os
import re

temp, density = "1e8", "1e3"

total_dir = f"../ppn/nuppn/frames/ppn/i_process_runs/t{temp}_p{density}"

frame_dir = f"ppn/nuppn/frames/ppn/i_process_runs/t{temp}_p{density}/ppn_frame.input"
abundance_dir = f"ppn/nuppn/frames/ppn/i_process_runs/t{temp}_p{density}/intershell_Z001 (1).txt"

def Temp_Density(filename):
    with open(filename, "r") as f:
        a = f.readlines()

    temp_line = a[2].split()[2]
    temperature = float(temp_line.split("d")[0]) * (10 ** (int(temp_line.split("d")[1]) + 9))

    density_line = a[3].split()[2]
    density = float(density_line.split("d")[0]) * (10 ** (int(density_line.split("d")[1])))
    return temperature, density


def fix_ele(line):
    elem = line[0:10].strip()
    reg = r"^(?P<elem>[a-zA-Z]+)\s*(?P<mass>\d+)"
    reg = re.compile(reg)
    e = re.fullmatch(reg, elem)
    if e is None:
        return None
    return e['elem'].capitalize() + " " + e['mass']

def abundance(filename, element):

    with open(filename, "r") as f:
        a = f.readlines()

    for num, line in enumerate(a):
        file_ele_name = fix_ele(a[num])
        # print(file_ele_name)
        if file_ele_name == element:
            return line.split()[2]


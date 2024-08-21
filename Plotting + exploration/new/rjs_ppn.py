# Analysis tools for ppn output

import numpy as np
import matplotlib.pyplot as plt
import re
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.patches import Arrow
from itertools import cycle

# Read in Solar abundance data
solar_file = open('../Data_files/iniab1.4E-02As09.ppn','r')

SolarData = []
PROT = solar_file.readline()
SolarData.append(PROT.split())
for line in solar_file.readlines()[1:]:
    # line = re.findall('\S+',line)
    p_num = line[:3].strip()
    name = line[3:6].strip()
    total = line[6:9].strip()
    n_num = int(total) - int(p_num)
    actual_name = name.strip()
    SolarData.append([p_num, name, total, line.split()[-1]])
solar_file.close()

# Run through Solar data and produce a list of Stable isotopes
# Also set up a dictionary for Solar Abundances
Stable_Isotopes = []
Solar_Abundance = {}
Solar_Element = {}
el_prev = "WTF"
X_tot = 0
for i in range(len(SolarData)):
# Sort out the inconsistent formatting of the data
    if SolarData[i][1] == "PROT":
        SolarData[i].append(SolarData[i][2])
        SolarData[i][2] = '1'
        SolarData[i][1] = "H"
    iso = SolarData[i][1].capitalize() + '-' + SolarData[i][2]
    el = SolarData[i][1].capitalize()
    Stable_Isotopes.append(iso)
    Solar_Abundance[iso] = float(SolarData[i][3])
    if el != el_prev:
        Solar_Element[el_prev] = X_tot
        el_prev = el
        X_tot = float(SolarData[i][3])
    else:
        X_tot = X_tot + float(SolarData[i][3])
# and don't forget to write the last element to the dictionary
Solar_Element[el] = X_tot

# create dictionary of element names from Z, Z from element names. Uses isotopesdatabase_cf.txt
Z_previous = 999
Element_Symbol = {}
Z_dict = {}

isotope_database = open('../ppn/nuppn/frames/ppn/i_process_runs/i_process_template/isotopedatabase_cf.txt', 'r')

while True:
# read in the isotope database, line by line
    line = isotope_database.readline()
    if not line:
        break
    line = re.findall('\S+',line)
    if len(line) == 0:
        line = "#"
# and only do something if it isn't the header part of the file
    if line[0] != "#":
        Z = line[0]
        el = line[2].capitalize()
        A = line[1]
        isotope = el + '-' + A
# Only add to the dictionaries if we found a new element
        if Z_previous != Z:
            Element_Symbol[int(Z)] = el
            Z_dict[el] = int(Z) 

# create a list of the stable z and n's from the stable isotope list in SolarData. Needed for plotting stable isotopes as solid squares in plot.
n_stable = []
z_stable = []
for i in range(len(SolarData)):
    iso = SolarData[i][1].capitalize()
    A = int(SolarData[i][2])
    n_stable.append(A - Z_dict[iso])
    z_stable.append(Z_dict[iso])

            
def read_xtime(file_path):
    # Read in the header of a file and neaten the isotope names up
    f = open(file_path,"r")

    headers = f.readline()
    headers = " ".join(headers.split())
    headers = headers.split("|")
    headers.pop(0)

    for i in range(len(headers)):
        if i == 0:
            headers[i] = "cycle"
        if i == 1:
            headers[i] = "time"
        if i == 2:
            headers[i] = "t9"
        if i == 3:
            headers[i] = "rho"
        if i == 4:
            headers[i] = "1-sum(yps)"
        if i == 5:
            headers[i] = "ye"
        if i == 6:
            headers[i] = "n"
        if i == 7:
            headers[i] = "H-1"
        if i > 7:
            el = re.findall("[A-Z]+",headers[i])[0]
            el = el.casefold()
            el = el.capitalize()
            A = re.findall("[0-9]+",headers[i])[-1]
            # currently isomers are labelled by letter in the 3rd position of the isomer name
            # but only with a lower case if it's less than the 26th level -- the below will break, eventually...
            level = re.findall("[a-z]",headers[i])
            if (len(level) == 0):
                iso = el + "-" + A
            else:
                iso = el + "-" + A + level[0]
            headers[i] = iso
    
    f.close()
    
    # now genfromtxt the actual data
    tmp_data = np.genfromtxt(file_path,skip_header=1)
    
    return tmp_data, headers
    
def reorder_isotopes(tmp_data,tmp_headers,tmp_cycle):
    # the damned output isn't in order, so let's reorder it...
    ZA = np.zeros(shape=(100,300))
    ZA_mask = np.zeros(shape=(100,300))

    for iso in tmp_headers:
        i = tmp_headers.index(iso)
        if iso !="cycle" and iso != 'time' and iso !='t9' and iso != "1-sum(yps)" and iso != "ye" and iso != "rho":
            if iso == "n":
                Z = 0
                A = 1
                ZA[Z][A] = tmp_data[tmp_cycle][i]
                ZA_mask[Z][A] = 1
            else:
                el = re.findall("[A-Za-z]+",iso)[0]
                A = int(re.findall("[0-9]+",iso)[0])
                Z = Z_dict[el]
                ZA[Z][A] = tmp_data[tmp_cycle][i]
                ZA_mask[Z][A] = 1
    # This is broken for isomers -- but luckily you get away with it: .index finds the FIRST instance of an occurrence, and the isomers are always at the end of the array. You really should fix this...

    new_model_headers = []
    new_model_headers.append("cycle")
    new_model_headers.append("time")
    new_model_headers.append("t9")
    new_model_headers.append("1-sum(yps)")
    new_model_headers.append("ye")
    new_model_headers.append("rho")

    reordered = []
    for i in range(0,6):
        reordered.append(tmp_data[tmp_cycle][i])
    for Z in range(0,100):
        for A in range(0,300):
            if ZA_mask[Z][A] == 1:
                Symbol=Element_Symbol.get(Z)
                iso_name = str(Symbol)+"-"+str(A)
                new_model_headers.append(iso_name)
                reordered.append(ZA[Z][A])
    return reordered, new_model_headers
    
def beta_decay_isotopes(tmp_data,tmp_headers,tmp_cycle):
    # take an input x-time cycle and decay its abundances to stable isotopes, based the Solar list
    # NB -- note this will royally screw up the heaviest isotopes where there are NO STABLE ELEMENTS!
    
    # create a stable isotope array in A,Z
    stable_mask = np.zeros(shape=(100,300))
    # populate with 1's for stable isotopes using n_stable, z_stable lists
    for n, Z in zip(n_stable, z_stable):
        A = n + Z
        stable_mask[Z][A] = 1
        
    X = np.zeros(shape=(100,300))
    ZA_mask = np.zeros(shape=(100,300))   
    # create an array of abundances from the x-time cycle
    for iso in tmp_headers:
        i = tmp_headers.index(iso)
        if iso !="cycle" and iso != 'time' and iso !='t9' and iso != "1-sum(yps)" and iso != "ye" and iso != "rho":
            if iso == "n":
                Z = 0
                A = 1
                X[Z][A] = tmp_data[tmp_cycle][i]
                ZA_mask[Z][A] = 1
            else:
                el = re.findall("[A-Za-z]+",iso)[0]
                A = int(re.findall("[0-9]+",iso)[0])
                Z = Z_dict[el]
                X[Z][A] = tmp_data[tmp_cycle][i]
                ZA_mask[Z][A] = 1
                
    # Go along each A, starting with the lowest Z, and move abundances to the next highest Z if the element is unstable
    
    for a in range(0,300): # this is too long -- is there a way to shorten this?
        # first pass through the stable list to find the highest Z that is stable
        z_max = 0
        for z in range(0,100):
            if stable_mask[z][a] == 1:
                z_max = z
        for z in range(0,100):
            if ZA_mask[z][a] == 1: # only do this if the isotope is in our network -- might be a problem if we decay to something that doesn't exist...
                if stable_mask[z][a] == 0 and X[z][a] > 1e-99 and z < z_max:
                    X[z+1][a] = X[z+1][a] + X[z][a] # unstable isotope, so add abundances to the previous one
                    X[z][a] = 0
                # now deal with the beta-plus decays -- note there will be some ambiguity for inbetweeners -- check with JK to decide which these are and what to do
                elif stable_mask[z][a] == 0 and z > z_max and X[z][a]:
                    X[z_max][a] = X[z_max][a] + X[z][a]
                    X[z][a] = 0
    
    # write the array to a new list for output
    new_model_headers = []
    new_model_headers.append("cycle")
    new_model_headers.append("time")
    new_model_headers.append("t9")
    new_model_headers.append("1-sum(yps)")
    new_model_headers.append("ye")
    new_model_headers.append("rho")

    reordered = []
    for i in range(0,6):
        reordered.append(tmp_data[tmp_cycle][i])
    for Z in range(0,100):
        for A in range(0,300):
            if ZA_mask[Z][A] == 1:
                Symbol=Element_Symbol.get(Z)
                iso_name = str(Symbol)+"-"+str(A)
                new_model_headers.append(iso_name)
                reordered.append(X[Z][A])
    return reordered, new_model_headers
    
def plot_abundance_time(isotopes,tmp_data,tmp_headers):
    
    # define cycle list for colours and linestyles
    lines = ["-","--",":","-."]
    colours = ['r','g','b','magenta','cyan','orange','black']
    linecycler = cycle(lines)
    colourcycler = cycle(colours)
    
    for iso in isotopes:
        iso_index = tmp_headers.index(iso)
        el = re.findall("[A-Za-z]+",iso)[0]
        A = re.findall("[0-9]+",iso)[0]
        iso_label = "$^{" + A + "}$" + el
        plt.plot([x[1]*60*60*24*365 for x in tmp_data],[x[iso_index] for x in tmp_data],label=iso_label,linestyle=next(linecycler),c=next(colourcycler))
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Time (s)")
    plt.ylabel("Mass fraction")
    
    plt.legend()
    
    return
    
def plot_isotopes(tmp_data,tmp_headers,tmp_cycle,**kwargs):
    # make an isotope plot, either as raw X, production factors (relative to Solar), or [A/B]
    
    plot_type = kwargs.get("plot",None)
    decay = kwargs.get("decay",True)
    if plot_type == "[A/B]":
        ref_isotope = kwargs.get("reference","Fe-56")
    
    if decay == True:
    # decay all isotopes
        X_cycle, new_headers = beta_decay_isotopes(tmp_data,tmp_headers,tmp_cycle)
    else:
        X_cycle, new_headers = reorder_isotopes(tmp_data,tmp_headers,tmp_cycle)
        
    xlim = kwargs.get("xlim",None)
    ylim = kwargs.get("ylim",None)
    
    fig, ax = plt.subplots()

    lines = ["-","--","-.",":"]
    linecycler = cycle(lines)
    
    markers = ["+","d","x","o","s"]
    markerscycler = cycle(markers)

    # set-up list before we process data
    A_array = []
    X_array = []
    el_prev = "WTF"
    
    
    # if using [A/B], pull the reference isotope and its initial value
    if plot_type == "[A/B]":
        ref_abundance = X_cycle[new_headers.index(ref_isotope)]
        ref_initial = Solar_Abundance[ref_isotope]
        ref_ratio = ref_abundance/ref_initial
    
    for i in range(len(new_headers)):
        if i > 6:
            iso = new_headers[i]
            X = X_cycle[i]
            el = re.findall("[a-zA-Z]+",iso)
            A = re.findall("\d+",iso)
            A = int(A[0])
        
        # are we on a new element? If so, plot the old pair, and reset the arrays
            if el != el_prev and len(A_array) > 0:
                plt.plot(A_array, X_array, marker=next(markerscycler), linestyle=next(linecycler))
                max_X = np.max(X_array)
                max_A = A_array[X_array.index(max_X)]
                ax.annotate(el_prev[0], (max_A,max_X),xytext=(0, 5),textcoords="offset points")
                A_array = []
                X_array = []
            
        
            if X > 0:
                A_array.append(A)
# process abundance data depending on type of plot requested
                if plot_type == "production":
        # check this is a stable element before trying to calculate things!
                    try:
                # get initial (Solar) abundance from Solar_Abundance dictionary
                        initial_abund = float(Solar_Abundance[iso])
                        X_array.append(X/initial_abund)
                    except:
                        X_array.append(0)
                elif plot_type == "[A/B]":
                    # check this is a stable element before trying to calculate things!
                    try:
                    # get initial (Solar) abundance from Solar_Abundance dictionary
                        initial_abund = float(Solar_Abundance[iso])
                        X_array.append(np.log10(X/initial_abund) - np.log10(ref_ratio))
                    except:
                        X_array.append(-99)
                else:
                    X_array.append(X)

            el_prev = el
            
    ax.grid(linestyle=":")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    plt.xlabel("Atomic Weight")

    # set labels depending on choice of plot
    if plot_type == "production":
        plt.ylabel("X/X$_\mathrm{i}$")
        plt.yscale('log')
    elif plot_type == "[A/B]":
        A = re.findall("\d+",ref_isotope)[0]
        el = re.findall("[a-zA-Z]+",ref_isotope)[0]
        y_label = "[X/$^{" + A + "}$" + el + "]"
        plt.ylabel(y_label)
    else:
        plt.ylabel("Mass fraction")
        plt.yscale('log')
    
    return
    
def plot_elements(tmp_data_list,tmp_headers_list,tmp_cycle_list,**kwargs):
    # NB -- to make this plot multiple models, ALL inputs must be in lists -- even if you only want to plot one data set, it still has to be in a list of one.
    
    # make an element plot, either as raw X, production factors (relative to Solar), or [A/B]
    # plot_abundance_z already does non-decayed elements, but could combine it here with a keyword 'decayed'
    plot_type = kwargs.get("plot",None)
    decay_in = kwargs.get("decay",[True])
    if plot_type == "[A/B]":
        ref_element = kwargs.get("reference","Fe")
    
    xlim = kwargs.get("xlim",None)
    ylim = kwargs.get("ylim",None)
    
    # set-up for plotting
    fig, ax = plt.subplots()

    lines = ["-","--","-.",":"]
    linecycler = cycle(lines)
    
    markers = ["+","d","x","o","s"]
    markerscycler = cycle(markers)
    
    # expecting an input list for decays, but maybe we didn't give one. If so, treat everything the same
    if len(decay_in) != len(tmp_cycle_list):
        i_one_decay = 1
    else:
        i_one_decay = 0
    
    for data_set in range(len(tmp_cycle_list)):
# ideally, decay should be in this list of lists as well, so we can plot decayed versus undecayed abundances...
        print(tmp_data_list[data_set])
        tmp_data = tmp_data_list[data_set]
        tmp_headers = tmp_headers_list[data_set]
        tmp_cycle = int(tmp_cycle_list[data_set]) # Oh python how I loathe your inability to covert types...
        if i_one_decay == 1:
            decay = decay_in[0]
        else:
            decay = decay_in[data_set]
        # Note I've assumed either all values or only one is given. If you've mismatched your input lists, you're on your own.
        
        if decay == True:
        # decay all isotopes
            X_cycle, new_headers = beta_decay_isotopes(tmp_data,tmp_headers,tmp_cycle)
            decay_label = ", decayed"
        else:
            X_cycle, new_headers = reorder_isotopes(tmp_data,tmp_headers,tmp_cycle)
            decay_label = ""
        
        # compress the isotope arrays and headers into elements
    
        X_element = []
        element_header = []
        Z = []
        el_prev = 'H'
#        el_prev = "NOT_AN_ELEMENT"
        X_tot = 0
        for i in range(len(new_headers)):
            if i > 6:
                iso = new_headers[i]
                X = X_cycle[i]
                el = re.findall("[a-zA-Z]+",iso)[0]
                if el != el_prev:
                    # we've changed element, so append the previous one and its total abundance to the list
                    X_element.append(X_tot)
                    element_header.append(el_prev)
                    Z.append(Z_dict[el_prev])       
                    el_prev = el
                    X_tot = X
                else:
                    X_tot = X_tot + X
        # ensure the last element gets added to the list
        X_element.append(X_tot)
        element_header.append(el)
        Z.append(Z_dict[el])
    
    # process the X array according to the type of plot we want
        if plot_type == "production":
            for i in range(len(X_element)):
                try:
                    X_element[i] = X_element[i]/Solar_Element[element_header[i]]
                except:
                    X_element[i] = 0
        elif plot_type == "[A/B]":
# work out np.log10 for X_ref/X_ref_i
            X_ref = X_element[element_header.index(ref_element)]
            X_ref_i = Solar_Element[ref_element]
            fac = np.log10(X_ref/X_ref_i)
            for i in range(len(X_element)):
                try:
                    X_element[i] = np.log10(X_element[i]/Solar_Element[element_header[i]]) - fac
                except:
                    X_element[i] = -99
    
        # set a label for the data
        label = str(tmp_cycle) + decay_label
        # now plot the data
        plt.plot(Z, X_element, marker=next(markerscycler), linestyle=next(linecycler), label=label)
                   
    ax.grid(linestyle=":")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    plt.xlabel("Proton Number")
    plt.legend()
    
    # set labels depending on choice of plot
    if plot_type == "production":
        plt.ylabel("X/X$_\mathrm{i}$")
        plt.yscale('log')
    elif plot_type == "[A/B]":
        el = re.findall("[a-zA-Z]+",ref_element)[0]
        y_label = "[X/" + el + "]"
        plt.ylabel(y_label)
    else:
        plt.ylabel("Mass fraction")
        plt.yscale('log')
    
    return
    
def plot_nuclide_chart(tmp_data,nmod,tmp_headers,**kwargs):
    xlim = kwargs.get("xlim",(0,10))
    ylim = kwargs.get("ylim",(0,10))
    clim = kwargs.get("clim",(-15,0))
    tmp_flux_data = kwargs.get("flux",[0])
    flux_scale = kwargs.get("scale","log")
    upper_flux = kwargs.get("upper_flux",None)
    lower_flux = kwargs.get("lower_flux",1e-20)
    
    # This routine assumes that Solar Data has been read in, as well as isotopedatabase.txt, in order to create the
    # isotope list for plotting get a list of n, z and abundance from the x-time file for all isotopes in the
    # network, for a chosen model number set model number here
    model_number = nmod # would be better to do with an index... or np.where, but its slow.

    n_list = []
    z_list = []
    abund = []
    for iso in tmp_headers:
        if iso !="cycle" and iso != 'time' and iso !='t9' and iso != "1-sum(yps)" and iso != "ye" and iso != "rho":
            if iso != "n":
                el = re.findall('[A-Za-z]+',iso)[0]
                A = int(re.findall('[0-9]+',iso)[0])
            else:
                el = "Nn"
                A = 1
            neutrons = A - Z_dict[el]
            i = tmp_headers.index(iso)
            abund.append(np.log10(tmp_data[model_number][i] + 1e-99))
            n_list.append(neutrons)
            z_list.append(Z_dict[el])
            

            # This is the actual plotting routine for the nuclide chart

    stablelist = []
    # create patches for stable isotopes
    for n, z in zip(n_stable, z_stable):
        rect = Rectangle((n,z),1,1,fill=None)
        stablelist.append(rect)
    
    stable = PatchCollection(stablelist,edgecolor="black",match_original=True,lw=3)
    
    # create patches for actual data

    patchlist = []
    for n, z in zip(n_list, z_list):
        rect = Rectangle((n,z),1,1)
        patchlist.append(rect)

    fig, ax = plt.subplots(figsize=(10,6))

    # currently this is 'plotting' all the nuclides -- we should have it only plot those in range
    for r in range(len(patchlist)):
        rx, ry = patchlist[r].get_xy()
        cx = rx + patchlist[r].get_width()/2.0
        cy = ry + patchlist[r].get_height()/2.0
    #    if rx != 0:
        if ry != 0:
            el = Element_Symbol.get(ry)
    #        iso_name = el +'-'+ str(int(rx+ry))
            iso_name = "$^{" + str(int(rx+ry)) + "}$" + el
        else:
                # we're dealing with neutrons that don't have a Z
            iso_name = "n"
        ax.annotate(iso_name, (cx, cy), color='black', weight='bold', 
                fontsize=12, ha='center', va='center')
    
    pc = PatchCollection(patchlist,edgecolor="black",cmap=plt.cm.hot_r)
    # set the abundance scale with clim
    pc.set_clim(clim)
    pc.set_array(np.array(abund))

    # choose which bits of the plot to show...
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('Neutrons')
    ax.set_ylabel('Protons')
    ax.add_collection(pc)
    ax.add_collection(stable)

    cbar = plt.colorbar(pc,orientation = "vertical")
    cbar.set_label("log X$_\mathrm{i}$")
    
    # plot flux arrows if flux data is present. Assumes a preset value of 1e-90 for flux cutoff
    if len(tmp_flux_data) != 1:
        # call routine to reduce flux data
        tmp_flux_arrows = fluxes_for_plotting(tmp_flux_data, 1e-50)
        # then call the flux arrow plotting routine
        # set a flux limit
        flux_limit = lower_flux #1e-20
        # and find the maximum flux (or set a max if you want to do it that way...)
        if upper_flux == None:
            max_flux = max([flux[4] for flux in tmp_flux_arrows])
        else:
            max_flux = 1e-15
        # scale for fluxes to be plotted on
        df = np.log10(max_flux/flux_limit)

        arrows = []
        flux_array = []
        for flux_data in tmp_flux_arrows:
            if abs(flux_data[4]) > flux_limit: # only plot it if its significant
                x = flux_data[1] + 0.5 # initial neutrons - with an offset to get box centre
                y = flux_data[0] + 0.5 # initial protons - with an offset to get box centre
                dx = flux_data[3] - flux_data[1] # difference in protons
                dy = flux_data[2] - flux_data[0]
#        print(x,y,dx,dy,flux_data[4])
                weight = np.log10(abs(flux_data[4])/flux_limit)/df
                arrow = Arrow(x, y, dx, dy, width=weight)
                arrows.append(arrow)
            # append flux to array for coloured arrow plotting!
                if flux_scale == "log":
                    flux_array.append(np.log10(abs(flux_data[4])+flux_limit/10))
                elif flux_scale == "linear":
                    flux_array.append(abs(flux_data[4]))
                else:
                    print("I don't know what you expect me to do...")
    
        flux_arrows = PatchCollection(arrows,edgecolor="black",cmap = plt.cm.jet)
        # set the abundance scale with clim
        flux_arrows.set_clim([np.log10(flux_limit), None])
        flux_arrows.set_array(np.array(flux_array))

        ax.add_collection(flux_arrows)
        
    return
    
def fluxes_for_plotting(temp_flux_data,flux_cutoff):
    # expects to be given an *array* of flux data from a flux_****.DAT file
    # returns an array that can be used to plot reaction arrows on a nuclide chart

    # first reduce the list to only relevant reactions, as defined by flux_cutoff
    reduced_flux = []
    for i in range(len(temp_flux_data)):
        if temp_flux_data[i][9] > flux_cutoff:
            tmp_reaction_data = []
            for j in range(0,10):
                tmp_reaction_data.append(temp_flux_data[i][j])
            reduced_flux.append(tmp_reaction_data)
    
    # reduce the flux list to net fluxes

    # find reaction pairs
    # first create a mask for the array
    paired = np.full(len(reduced_flux),False)

    pair_list = []
    # go through the array and look for pairs
    for i in range(len(reduced_flux)):
        # Only check if we haven't paired this reaction yet
        if paired[i] == False:
            for j in range(len(reduced_flux)):
            # Only check if we haven't paired this reaction yet
                if paired[j] == False:
                # do the Z1,A1,Z2,A2 match the Z3,A3,Z4,A4
                # and I can only hope they've kept the light - heavy order in both...
                # aannnddd they haven't so we need extra logic...
                    if (reduced_flux[i][1] == reduced_flux[j][5] and reduced_flux[i][2] == reduced_flux[j][6] and reduced_flux[i][3] == reduced_flux[j][7] and reduced_flux[i][4] == reduced_flux[j][8]) or (reduced_flux[i][1] == reduced_flux[j][7] and reduced_flux[i][2] == reduced_flux[j][8] and reduced_flux[i][3] == reduced_flux[j][5] and reduced_flux[i][4] == reduced_flux[j][6]):     
                        # and the inverse pairing also checking both possible orders
                        if (reduced_flux[j][1] == reduced_flux[i][5] and reduced_flux[j][2] == reduced_flux[i][6] and reduced_flux[j][3] == reduced_flux[i][7] and reduced_flux[j][4] == reduced_flux[i][8]) or (reduced_flux[j][1] == reduced_flux[i][7] and reduced_flux[j][2] == reduced_flux[i][8] and reduced_flux[j][3] == reduced_flux[i][5] and reduced_flux[j][4] == reduced_flux[i][6]):   
                            #print("Pair found!", i, j)
                            #print(reduced_flux[i][1:5],reduced_flux[i][5:9])
                            #print(reduced_flux[j][5:9],reduced_flux[j][1:5])
                            pair = []
                            pair.append(i)
                            pair.append(j)
                            pair_list.append(pair)
                            paired[i] = True
                            paired[j] = True

    # Use the paired list to calculate net flux
    flux_plot_data = []
    for pair in pair_list:
        i = pair[0]
        j = pair[1]
        net_flux = reduced_flux[i][9] - reduced_flux[j][9]
        if net_flux < 0:
            # write to an array the pairing j-values, i-values and flux
            # I only want the Z,N of the heavier species
            flux_details = []
            flux_details.append(reduced_flux[j][1])
            flux_details.append(reduced_flux[j][2]-reduced_flux[j][1])
            flux_details.append(reduced_flux[i][1])
            flux_details.append(reduced_flux[i][2]-reduced_flux[i][1])
            flux_details.append(net_flux)
            flux_plot_data.append(flux_details)
        elif net_flux > 0:
            flux_details = []
            flux_details.append(reduced_flux[i][1])
            flux_details.append(reduced_flux[i][2]-reduced_flux[i][1])
            flux_details.append(reduced_flux[j][1])
            flux_details.append(reduced_flux[j][2]-reduced_flux[j][1])
            flux_details.append(net_flux)
            flux_plot_data.append(flux_details)
            # write to an array the pairing i-values, j-values and flux
        # and if it was zero, just don't both with it.
        
    # now use the paired list to find all the unpaired data
    for i in range(len(reduced_flux)):
        if paired[i] == False and reduced_flux[i][9] != 1e-99:
            # check if the flux is positive or negative -- and ignore it if its 1e-99 which is effectively zero
            # "Zero is less than nothing" -- M. Pignatari
            if reduced_flux[i][9] > 0:
                # positive flux, so we can take the Z, A values from the first column and the last two columns
                flux_details = []
                flux_details.append(reduced_flux[i][1])
                flux_details.append(reduced_flux[i][2]-reduced_flux[i][1])
                flux_details.append(reduced_flux[i][7])
                flux_details.append(reduced_flux[i][8]-reduced_flux[i][7])
                flux_details.append(reduced_flux[i][9])
                flux_plot_data.append(flux_details)
            else:
                # negative flux, so we can take the Z, A values from the last two column and the first two columns
                flux_details = []
                flux_details.append(reduced_flux[i][7])
                flux_details.append(reduced_flux[i][8]-reduced_flux[i][7])
                flux_details.append(reduced_flux[i][1])
                flux_details.append(reduced_flux[i][2]-reduced_flux[i][1])
                flux_details.append(reduced_flux[i][9])
                flux_plot_data.append(flux_details)

    return flux_plot_data
    

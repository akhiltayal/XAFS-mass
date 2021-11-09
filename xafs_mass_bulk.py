# -*- coding: utf-8 -*-
"""
Created on Sat Oct  2 14:18:07 2021

@author: tayalakh
"""

# %matplotlib notebook
# import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd
path_to_xafs_mass = "D:/Software/XAFSmass-1.4.0/kklmn-XAFSmass-fc099e2/XAFSmass/"
sys.path.append(path_to_xafs_mass)
# sys.path("D:/Software/XAFSmass-1.4.0/kklmn-XAFSmass-fc099e2/XAFSmass/")
from XAFSmassCalc import *

# from IPython.display import display

def powder(samples, mud=2, die_diam=10, exclude_elements = None, edge_to_find = None):
    """Xafs mass calculator for multiple samples"""
    exclu_elem = ["H", "He", "Li", "Be", "B", "C", "O", "N", "F",
             "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca" ]
    if exclude_elements is None:
        pass
    else:
        exclu_elem.extend(exclude_elements)  
    b = []
    for sample in samples:
        res = formula.parseString(sample).asList()
        for i in range(len(res)):
            d = []
            if res[i][0] in exclu_elem:
                continue
            else:
                c = res[i][0]+str(res[i][1])
                for j in range(len(res)):
                    if j==i:
                        continue
                    else:
                        c +=res[j][0]+str(res[j][1])
            d.append(c)
            with open(path_to_xafs_mass+'/data/Energies.txt', 'r') as en:
                energy = en.readlines()
                for line in energy:
                    elm_info = line.split()
                    if elm_info[1] == res[i][0] and edge_to_find is None and int(elm_info[0])>=50:
                        d.append(int(float(elm_info[5]))+50)
                    elif elm_info[1] == res[i][0] and edge_to_find is None and int(elm_info[0])<=57:
                        d.append(int(float(elm_info[2]))+50)
                    elif elm_info[1] == res[i][0] and edge_to_find=='K':
                        d.append(int(float(elm_info[2]))+50)
                    elif elm_info[1] == res[i][0] and edge_to_find=='L3':
                        d.append(int(float(elm_info[5]))+50)
            b.append(d) 
    tests = b
    column_names = ['Sample', 'Edge[eV]', 'E+50[eV]', 'µd', 'pellet_dia[mm]', 'nu[mmol]', 'mass[mg]', 'thickness[µm]', 
                'jump', 'Cellulose[mg]', 'BN[mg]']
    df = pd.DataFrame(columns = column_names)
    i = 0
    area = (die_diam*0.1/2)**2*np.pi
    for comp, E in tests:
        results = formula.parseString(comp)
        nu, m, th, eDict = calculate_powder(results.asList(), E, mud, area,
                                            rho=2.1)
        eDict2 = dict([k, [v[3], v[-2], v[-1]]] for k, v in eDict.items())
        massCell = area*0.1*1500
        massBN = area*0.1*2100
#         print(u'{0} at {1}eV: nu={2}mmol, mass={3}mg, thickness={4}µm, jump={5}'.format(
#               comp, E, round_to_n(nu), round_to_n(m), round_to_n(th), round_to_n(eDict2[results[0][0]][-1])))
        data = [comp, E-50, E, mud, die_diam, round_to_n(nu), round_to_n(m), round_to_n(th), 
                round_to_n(eDict2[results[0][0]][-1]), round_to_n(massCell), round_to_n(massBN)] 
        df.loc[i] = data
        i +=1
# powder(b, mud=1.9, die_diam=8)
    new = df.sort_values('E+50[eV]')
#     indexnames = new[(new.Sample.str[0]=='O') | (new.Sample.str[0]=='N') | (new.Sample.str[0]=='H') | (new.Sample.str[0]=='S') | (new.Sample.str[0]=='C')  ].index

#     indexnames = new[(new.Sample.str[0]=='H')].index
    
    
    # indexnames = new[(new.Sample.str[0]=='O') | (new.Sample.str[0]=='N') 
    #                  |(new.Sample.str[0:2]=='Mn') |(new.Sample.str[0:2]=='Fe') 
    #                  |(new.Sample.str[0:2]=='Ni') |(new.Sample.str[0:2]=='Co')
    #                 |(new.Sample.str[0:2]=='n')].index
    return(new)


    




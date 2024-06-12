# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 23:31:06 2020

@author: ashraya
"""

import numpy as np
import math
phi=[]
psi=[]
infile=open("/home/ashraya/MD_Project/my_md_files/jak2_human/dfg_relaxed/trial_3/4iva_relaxed_3_dfg_asp_phipsi.csv","r")
for line in infile:
    lineparts=line[0:-1].split(",")
    phi.append(float(lineparts[0]))
    psi.append(float(lineparts[1]))    
infile.close()
phi=np.radians(np.array(phi))
psi=np.radians(np.array(psi))
r_square=(np.square(np.sum(np.cos(phi)))+np.square(np.sum(np.sin(phi)))+np.square(np.sum(np.cos(psi)))+np.square(np.sum(np.sin(psi))))/2.0
r_av=math.sqrt(r_square)/len(phi)
circvar=1-r_av
print circvar
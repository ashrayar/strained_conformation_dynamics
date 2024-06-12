# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 00:53:28 2020

@author: ashraya
"""
from collections import Counter
infile=open("/home/ashraya/MD_Project/all_protein_nma/nma_loop_flucts_all_norm.csv","r")
strained_res=[]
relaxed_res=[]
strained_res_nogly=[]
relaxed_res_nogly=[]
for line in infile:
    lineparts=line[0:-1].split(",")
    if lineparts[-1] in ("Allowed","OUTLIER"):
        strained_res.append(lineparts[3])
        if "GLY" not in line:
            strained_res_nogly.append(lineparts[3])
    else:
        relaxed_res.append(lineparts[3])
        if "GLY" not in line:
            relaxed_res_nogly.append(lineparts[3])
        
infile.close()
print Counter(strained_res).keys()
print Counter(strained_res).values()
strained_res_dist=dict(Counter(strained_res))
relaxed_res_dist=dict(Counter(relaxed_res))
#strained_res_dist={'TYR': 299, 'ILE': 305, 'GLN': 167, 'SER': 455, 'LYS': 393, 'GLY': 602, 'ALA': 417, 'GLU': 315, 'CYS': 107, 'ASP': 975, 'PRO': 337, 'MET': 80, 'THR': 307, 'VAL': 525, 'ASN': 829, 'ARG': 347, 'PHE': 297, 'HIS': 264, 'TRP': 112, 'LEU': 228}
strained_res_dist={'TYR': 299, 'ILE': 305, 'GLN': 167, 'SER': 455, 'LYS': 393, 'ALA': 417, 'GLU': 315, 'CYS': 107, 'ASP': 975, 'PRO': 337, 'MET': 80, 'THR': 307, 'VAL': 525, 'ASN': 829, 'ARG': 347, 'PHE': 297, 'HIS': 264, 'TRP': 112, 'LEU': 228}
#relaxed_res_dist={'TYR': 4294, 'ILE': 4763, 'GLN': 5931, 'SER': 11990, 'LYS': 8902, 'GLY': 23067, 'ALA': 11828, 'GLU': 8201, 'CYS': 1742, 'HIS': 4135, 'PRO': 14218, 'ASP': 13470, 'MET': 2056, 'THR': 10935, 'VAL': 6506, 'ASN': 10359, 'ARG': 6479, 'PHE': 4758, 'TRP': 1905, 'LEU': 9253}
relaxed_res_dist={'TYR': 4294, 'ILE': 4763, 'GLN': 5931, 'SER': 11990, 'LYS': 8902, 'ALA': 11828, 'GLU': 8201, 'CYS': 1742, 'HIS': 4135, 'PRO': 14218, 'ASP': 13470, 'MET': 2056, 'THR': 10935, 'VAL': 6506, 'ASN': 10359, 'ARG': 6479, 'PHE': 4758, 'TRP': 1905, 'LEU': 9253}
strained_res_percent=dict()
for i in strained_res_dist:
    strained_res_percent[i]=strained_res_dist[i]*100/float(len(strained_res_nogly))
relaxed_res_percent=dict()
for i in relaxed_res_dist:
    relaxed_res_percent[i]=relaxed_res_dist[i]*100/float(len(relaxed_res_nogly))
import matplotlib.pyplot as plt
import numpy as np
relaxed_res_list=np.arange(len(relaxed_res_dist))
plt.bar(relaxed_res_list, strained_res_percent.values(), width=0.3,align='center',color='red',alpha=0.7)        
plt.bar(relaxed_res_list-0.3, relaxed_res_percent.values(), width=0.3,align='center',color='green',alpha=0.7)        
plt.xticks(range(len(strained_res_percent)), strained_res_percent.keys(),fontsize=16)
plt.ylim(0,16)
plt.yticks(fontsize=16)
plt.xlim(-1,19)
plt.grid(b=True,which='major',axis='y',linestyle=":")
plt.xlabel("Amino acid residues",fontsize=22)
plt.ylabel("Percent",fontsize=22)
plt.legend(("Strained Residues","Relaxed Residues"),fontsize=22,loc="upper left")
plt.savefig("/home/ashraya/MD_Project/plots/strained_relaxed_res_aa_distribution_new.png",dpi=300,bbox_inches="tight")

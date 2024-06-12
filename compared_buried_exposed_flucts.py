# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 02:16:48 2020

@author: ashraya
"""
import numpy as np
from scipy import stats
infile=open("/home/ashraya/MD_Project/all_protein_nma/max_flucts_raw_chains_termini_filtered.csv","r")
exposed_flucts=[]
buried_flucts=[]
for line in infile:
    if "nan\n" in line:
        continue
    flucts=dict()
    print line[0:4]
    fluctfile=open("/home/ashraya/MD_Project/all_protein_nma/fluctuations_raw/normalized_flucts/"+line[0:4]+"_norm_flucts.csv","r")
    for line1 in fluctfile:
        line1parts=line1.split(",")
        resno,chain,na=line1parts[0].split("_")
        flucts[resno+"_"+chain]=line1parts[1]
    fluctfile.close()
    sasafile=open("/home/ashraya/MD_Project/all_protein_nma/all_protein_sasa/"+line[0:4]+"_fixed.rsa","r")
    for line1 in sasafile:
        if line1.startswith("RES"):
            line1parts=line1.split()
            res_chain=line1parts[3]+"_"+line1parts[2]
            if res_chain in flucts.keys():
                sasa=float(line1parts[5])
                if sasa > 7.0:
                    exposed_flucts.append(float(flucts[res_chain]))
                else:
                    buried_flucts.append(float(flucts[res_chain]))
    sasafile.close()
infile.close()
print np.mean(exposed_flucts)
print np.mean(buried_flucts)
import matplotlib.pyplot as plt
ex_bins=np.linspace(-3,18,43)
weights = np.ones_like(exposed_flucts)/float(len(exposed_flucts))
x1,b1,p1=plt.hist(exposed_flucts,bins=ex_bins,weights=weights,alpha=0.5,label="exposed residues",color="green")
plt.axvline(np.mean(exposed_flucts),color="green",linestyle="dashed",linewidth=2)
weights = np.ones_like(buried_flucts)/float(len(buried_flucts))
bu_bins=np.linspace(-2,5,15)
x2,b2,p2=plt.hist(buried_flucts,weights=weights,bins=bu_bins,alpha=0.5,label="buried residues",color="red")
plt.axvline(np.mean(buried_flucts),color="red",linestyle="dashed",linewidth=2)
plt.legend(loc="upper right",fontsize=22)
#plt.xlim(-5,22)
#plt.xlim(-2,18)
plt.xticks(np.linspace(-2,18,21),fontsize=16)
plt.xlim(-2,13)
#plt.xticks(np.linspace(0,1,11),fontsize=16)
plt.yticks(fontsize=16)
plt.grid(b=True,which='major',axis='y',linestyle=":")
plt.xlabel("Normalized Square Fluctuations",fontsize=22)
plt.ylabel("Proportion of residues",fontsize=22)

combined=[exposed_flucts,buried_flucts]
fig,axis=plt.subplots()
bplot=axis.boxplot(combined,labels=["Exposed Residues","Buried Residues"],patch_artist=True,showmeans=True)
colors=['pink','lightgreen']
axis.set_xticklabels(["Exposed Residues","Buried Residues"],fontsize=22)
axis.tick_params(axis='y',which='major', labelsize=16)
for element in ['boxes','whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(bplot[element], color='black')  
for patch,color in zip(bplot['boxes'],colors):
    patch.set_facecolor(color)
axis.yaxis.grid(True)
axis.set_ylim(-3,18)
axis.set_ylabel("Normalized Fluctuations", fontsize=22)
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 16:24:48 2020

@author: ashraya
"""
import scipy.stats as stat
infile=open("/home/ashraya/MD_Project/my_md_files/aurka_human/5orl_relaxed/frustratometer_results/FrustrationData/5orl_relaxed_frustration.pdb_singleresidue","r")
x=infile.readline()
frust=[]
for line in infile:
    lineparts=line.split()
#    if lineparts[0] in ['159','160']:
#        continue
    frust.append(float(lineparts[-1]))
infile.close()
infile=open("/home/ashraya/MD_Project/my_md_files/aurka_human/5orl_relaxed/trial_3/trial3_rmsf.xvg","r")
rmsf=[]
for line in infile:
    if line.startswith(" "):
        lineparts=line.split()
        rmsf.append(float(lineparts[1]))
infile.close()
z_rmsf=stat.zscore(rmsf)

hf_rmsf=[]
nf_rmsf=[]
mf_rmsf=[]
for i in range(0,len(frust)):
    if frust[i]>=0.78:
        mf_rmsf.append(z_rmsf[i])
    elif frust[i]<-1:
        hf_rmsf.append(z_rmsf[i])
    else:
        nf_rmsf.append(z_rmsf[i])
print np.mean(hf_rmsf)
print np.mean(mf_rmsf)
print np.mean(nf_rmsf)
print stat.mannwhitneyu(hf_rmsf,mf_rmsf)
print stat.mannwhitneyu(hf_rmsf,nf_rmsf)
print stat.mannwhitneyu(mf_rmsf,nf_rmsf)
#print stat.pearsonr(z_rmsf,frust)
#print stat.spearmanr(z_rmsf,frust)

#import matplotlib.pyplot as plt
#import numpy as np
#plt.plot(z_rmsf)
#plt.plot(frust)
#plt.scatter(z_rmsf,frust)
#m, b = np.polyfit(z_rmsf, frust, 1)

#import seaborn as sns
#ax=sns.regplot(x=np.array(z_rmsf),y=np.array(frust))
#ax.set_xlabel("Normalized RMSF (from MD)",fontsize=22)
#ax.set_xticklabels(ax.get_xticks(),fontsize=16)
#ax.set_yticklabels(ax.get_yticks(),fontsize=16)
#ax.set_ylabel("Frustration Index",fontsize=22)
#plt.savefig("/home/ashraya/MD_Project/my_md_files/jak2_relaxed_trial3_frust.png",dpi=300,bbox_inches="tight")
        

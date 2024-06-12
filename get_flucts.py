# -*- coding: utf-8 -*-
"""
Created on Sat Jun  9 00:04:46 2018

@author: ashraya
"""
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats

infile=open("/home/ashraya/MD_Project/lower_res_analysis/all_kinases_dfg_molprob.csv","r")
unstrained_structs=[]
for line in infile:
    lineparts=line.split(",")
    if lineparts[-1] in ["Favored\n"]:#["Allowed\n","OUTLIER\n"]:
        unstrained_structs.append(lineparts[0])
infile.close()

infile=open("/home/ashraya/MD_Project/lower_res_analysis/nma_pdb_list_final.csv","r")
nma_structs=[]
for line in infile:
    nma_structs.append(line[0:-1])
infile.close()

#
#infile=open("/home/ashraya/MD_Project/Fhelix_asp_locations.csv","r")
#opfile=open("/home/ashraya/MD_Project/Fhelix-D_unstrained_flucts.csv","w")
#hrd_flucts=[]
#for line in infile:
#    lineparts=line.split(",")
#    if lineparts[0] in unstrained_structs:
#        residue=lineparts[1][0:-1]
#        infile1=open("fluctuations/"+lineparts[0]+"_unstrained_flucts.csv","r")
#        for line1 in infile1:
#            resid=line1.split(",")[0].split("_")[0]
#            if resid==residue:
#                hrd_fluct=(round(float(line1.split(",")[1][0:-1]),4))
#                opfile.write(lineparts[0]+","+resid+","+str(hrd_fluct)+"\n")
#                break
#        infile1.close()
#opfile.close()
#infile.close()
infile=open("/home/ashraya/MD_Project/lower_res_analysis/all_kinases_dfg_molprob.csv","r")
opfile=open("lower_res_analysis/all_kinases_dfg_unstrained_flucts.csv","w")
for line in infile:
    lineparts=line.split(",")
    residue=lineparts[1]
    if lineparts[0]  in unstrained_structs and lineparts[0] in nma_structs:
        infile1=open("lower_res_analysis/fluctuations/"+lineparts[0]+"_flucts.csv","r")
        for line1 in infile1:
            resid=line1.split(",")[0].split("_")[0]
            if resid==residue:
                hrd_fluct=(round(float(line1.split(",")[1][0:-1]),4))
                opfile.write(lineparts[0]+","+resid+","+str(hrd_fluct)+"\n")
                break
        infile1.close()
infile.close()
opfile.close()

#df=pd.read_csv("Fhelix-D_strained_flucts.csv",names=["structure","resno","fluctuation"])
df=pd.read_csv("/home/ashraya/MD_Project/lower_res_analysis/all_kinases_dfg_strained_flucts.csv",names=["structure","resno","fluctuation"])
strained_flucts=df['fluctuation']
#df=pd.read_csv("Fhelix-D_unstrained_flucts.csv",names=["structure","resno","fluctuation"])
df=pd.read_csv("/home/ashraya/MD_Project/lower_res_analysis/all_kinases_dfg_unstrained_flucts.csv",names=["structure","resno","fluctuation"])
unstrained_flucts=df['fluctuation']
combined=[strained_flucts,unstrained_flucts]

fig,axis=plt.subplots()
bplot=axis.boxplot(combined,labels=["DFG_Asp-Strained","DFG_Asp-Relaxed"],patch_artist=True)
#bplot=axis.boxplot(combined,labels=["Inactive","active"],patch_artist=True)

colors=['pink','lightgreen']
for element in ['boxes','whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(bplot[element], color='black')  
for patch,color in zip(bplot['boxes'],colors):
    patch.set_facecolor(color)
axis.yaxis.grid(True)
axis.set_ylabel("Normalized Fluctuations")
#stats.ttest_ind(strained_flucts,unstrained_flucts,equal_var=False)
stats.mannwhitneyu(strained_flucts,unstrained_flucts)


# New Normalization
infile=open("/home/ashraya/MD_Project/hrd_flucts_no_ligands_unstrained.csv","r")
opfile=open("/home/ashraya/MD_Project/hrd_flucts_unstrained_newnorm.csv","w")
for line1 in infile:
    line1parts=line1[0:-1].split(",")
    resno=line1parts[1]
    fluctfile=open("/home/ashraya/MD_Project/fluctuations_no_ligands/normalized_flucts/"+line1parts[0]+"_norm_flucts.csv","r")
    for line in fluctfile:
        resdetails,znorm,maxnorm=line[0:-1].split(",")
        res=resdetails.split("_")[0]
        if res==resno:
            opfile.write(line1parts[0]+","+res+","+znorm+","+maxnorm+"\n")
            break
    fluctfile.close()
infile.close()
opfile.close()

df=pd.read_csv("/home/ashraya/MD_Project/dfg_flucts_strained_newnorm.csv",names=["structure","resno","fluct_znorm","fluct_maxnorm"])
strained_flucts=df['fluct_znorm']
strained_flucts=df['fluct_maxnorm']
df=pd.read_csv("/home/ashraya/MD_Project/dfg_flucts_unstrained_newnorm.csv",names=["structure","resno","fluct_znorm","fluct_maxnorm"])
unstrained_flucts=df['fluct_znorm']
unstrained_flucts=df['fluct_maxnorm']
print np.mean(strained_flucts)
print np.mean(unstrained_flucts)
stats.mannwhitneyu(strained_flucts,unstrained_flucts)
combined=[strained_flucts,unstrained_flucts]
fig,axis=plt.subplots()
bplot=axis.boxplot(combined,labels=["DFG_Asp-Strained","DFG_Asp-Relaxed"],patch_artist=True,showmeans=True)
#bplot=axis.boxplot(combined,labels=["HRD_Arg-Strained","HRD_Arg-Relaxed"],patch_artist=True,showmeans=True)
colors=['pink','lightgreen']
axis.tick_params(axis='y',which='major', labelsize=16)
axis.set_xticklabels(["Strained Residues","Relaxed Residues"],fontsize=22)
for element in ['boxes','whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(bplot[element], color='black')  
for patch,color in zip(bplot['boxes'],colors):
    patch.set_facecolor(color)
axis.yaxis.grid(True)
axis.set_ylabel("Normalized Fluctuations",fontsize=22)
#stats.ttest_ind(strained_flucts,unstrained_flucts,equal_var=False)
stats.mannwhitneyu(strained_flucts,unstrained_flucts)          
        
    

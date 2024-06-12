# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 23:55:39 2019

@author: ashraya
"""
'''
import MySQLdb
db = MySQLdb.connect("localhost","root","admin","KINASE_PROJECT" )
cursor = db.cursor() 

infile=open("/home/ashraya/MD_Project/multiple_occ_analysis/structures_1.0_june_2019.txt","r")
#infile=["5wqj\n"]
opfile=open("/home/ashraya/MD_Project/loop_residue_bfactors_1A.csv","w")
for line1 in infile:
    if line1.startswith("#"):
        continue
    print line1[0:-1]
    dsspfile=open("/home/ashraya/MD_Project/dssp_outputs/"+line1[0:-1]+".dssp","r")
    loops=[]
    flag=0
    for line in dsspfile:
        if line.startswith("  #  RESIDUE"):
            flag=1
            continue
        if flag==1:
            if "!*" in line:
                continue
            dssp=line[16]
            if dssp in ("S","T"," "):
                resno=line[5:10].strip()
                chain=line[11]
                loops.append(resno+","+chain)
    dsspfile.close()
    query="select chain_id,resno,resname,bfactor,molprob from protein_strain_bfactors where pdb_id='"+line1[0:-1]+"';"    
    cursor.execute(query)
    results=cursor.fetchall()
    for item in results:
        resno=str(item[1])
        to_search=resno+","+item[0]
        if to_search in loops:
            opfile.write(line1[0:-1]+","+item[0]+","+resno+","+item[2]+","+str(item[3])+","+item[4]+"\n")
opfile.close()
infile.close()'''


#Identifying alpha helices and beta sheet residues  
import MySQLdb
db = MySQLdb.connect("localhost","root","admin","KINASE_PROJECT" )
cursor = db.cursor() 

infile=open("/home/ashraya/MD_Project/multiple_occ_analysis/structures_1.0_june_2019.txt","r")
#infile=["5wqj\n"]
opfile=open("/home/ashraya/MD_Project/loop_residue_bfactors_1A.csv","w")
for line1 in infile:
    if line1.startswith("#"):
        continue
    print line1[0:-1]
    dsspfile=open("/home/ashraya/MD_Project/dssp_outputs/"+line1[0:-1]+".dssp","r")
    loops=[]
    flag=0
    for line in dsspfile:
        if line.startswith("  #  RESIDUE"):
            flag=1
            continue
        if flag==1:
            if "!*" in line:
                continue
            dssp=line[16]
            if dssp in ("S","T"," "):
                resno=line[5:10].strip()
                chain=line[11]
                loops.append(resno+","+chain)
    dsspfile.close()
    query="select chain_id,resno,resname,bfactor,molprob from protein_strain_bfactors where pdb_id='"+line1[0:-1]+"';"    
    cursor.execute(query)
    results=cursor.fetchall()
    for item in results:
        resno=str(item[1])
        to_search=resno+","+item[0]
        if to_search in loops:
            opfile.write(line1[0:-1]+","+item[0]+","+resno+","+item[2]+","+str(item[3])+","+item[4]+"\n")
opfile.close()
infile.close()
import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt

'''
loop_df=pd.read_csv("/home/ashraya/MD_Project/loop_residue_bfactors.csv",names=["pdb_id","chain","resno","resname","bfactor","molprob"])

outlier_df_1a=loop_df.loc[loop_df["molprob"].isin(["OUTLIER","Allowed"])]
relax_df_1a=loop_df.loc[loop_df["molprob"]=="Favored"]
outlier_bfac_1a=outlier_df_1a['bfactor']
relax_bfactor_1a=relax_df_1a['bfactor']
np.mean(outlier_bfac_1a)
np.mean(relax_bfactor_1a)
combined_1a=[outlier_bfac_1a,relax_bfactor_1a]
fig,axis=plt.subplots()
bplot=axis.boxplot(combined_1a,labels=["Strained Residues","Relaxed Residues"],patch_artist=True,showmeans=True)
colors=['pink','lightgreen']
axis.set_xticklabels(["Strained Residues","Relaxed Residues"],fontsize=22)
axis.tick_params(axis='y',which='major', labelsize=16)
for element in ['boxes','whiskers', 'fliers', 'means', 'medians', 'caps']:
    plt.setp(bplot[element], color='black')  
for patch,color in zip(bplot['boxes'],colors):
    patch.set_facecolor(color)
axis.yaxis.grid(True)
axis.set_ylabel("Normalized B-Factors", fontsize=22)
bins=range(-4,22,1)
x1,b1,p1=plt.hist(relax_bfactor_1a,bins=bins,normed=True,alpha=0.5,label="relaxed residues",color="green")
plt.axvline(np.median(relax_bfactor_1a),color="green",linestyle="dashed",linewidth=2)
x2,b2,p2=plt.hist(outlier_bfac_1a,bins=bins,normed=True,alpha=0.5,label="strained residues",color="red")
plt.axvline(np.median(outlier_bfac_1a),color="red",linestyle="dashed",linewidth=2)
plt.legend(loc="upper right")
plt.xlim(-5,22)
plt.xticks(range(-5,22,1),fontsize=14)
plt.yticks(fontsize=14)
plt.grid(b=True,which='major',axis='y',linestyle=":")
plt.xlabel("Normalized B-Factor (Strained Residues)",fontsize=20)
plt.ylabel("Normalized Frequency",fontsize=20)


plt.hist([outlier_bfac_1a,relax_bfactor_1a],bins=bins,normed=True)
plt.legend(loc="upper right")'''


strain_flucts_maxnorm=[]
unstrain_flucts_maxnorm=[]
strain_flucts_znorm=[]
unstrain_flucts_znorm=[]
#infile=open("/home/ashraya/MD_Project/all_protein_nma/nma_loop_flucts_1.2.csv","r")
infile=open("/home/ashraya/MD_Project/all_protein_nma/nma_loop_flucts_all_norm.csv","r")
for line in infile:
    lineparts=line.split(",")
    if lineparts[6][0:-1] in ("OUTLIER","Allowed"):
        strain_flucts_maxnorm.append(float(lineparts[5]))
        strain_flucts_znorm.append(float(lineparts[4]))
    elif lineparts[6][0:-1] in ("Favored"):
        unstrain_flucts_maxnorm.append(float(lineparts[5]))
        unstrain_flucts_znorm.append(float(lineparts[4]))
infile.close()

#import random
#opfile=open("/home/ashraya/MD_Project/all_protein_nma/sampling_comparison_outlier_favored_1.2.csv","w")
#for i in range(0,437):
#    unstrain_sample=random.sample(unstrain_flucts,424)
#    opfile.write(str(i+1)+','+str(stats.ttest_ind(strain_flucts,unstrain_sample,equal_var=False)[1])+'\n')
#opfile.close()
    

combined=[strain_flucts_znorm,unstrain_flucts_znorm]
#bins=np.linspace(0,1,21)
bins=np.linspace(-2,18,41)
#bins=np.linspace(-2,13,31)
weights = np.ones_like(unstrain_flucts_znorm)/float(len(unstrain_flucts_znorm))
x1,b1,p1=plt.hist(unstrain_flucts_znorm,bins=bins,weights=weights,alpha=0.5,label="relaxed residues",color="green")
plt.axvline(np.mean(unstrain_flucts_znorm),color="green",linestyle="dashed",linewidth=2)
weights = np.ones_like(strain_flucts_znorm)/float(len(strain_flucts_maxnorm))
x2,b2,p2=plt.hist(strain_flucts_znorm,weights=weights,bins=bins,alpha=0.5,label="strained residues",color="red")
plt.axvline(np.mean(strain_flucts_znorm),color="red",linestyle="dashed",linewidth=2)
plt.legend(loc="upper right",fontsize=22)
#plt.xlim(-5,22)
plt.xlim(-2,18)
#plt.xticks(np.linspace(-2,18,21),fontsize=16)
plt.xticks(np.linspace(-2,18,11),fontsize=16)
plt.xlim(-2,13)
plt.xticks(np.linspace(-2,14,9),fontsize=16)
#plt.xticks(np.linspace(0,1,11),fontsize=16)
plt.yticks(fontsize=16)
plt.grid(b=True,which='major',axis='y',linestyle=":")
plt.xlabel("Normalized Square Fluctuations",fontsize=22)
plt.ylabel("Proportion of residues",fontsize=22)
#fig,axis=plt.subplots()
#bplot=axis.boxplot(combined,labels=["Strained Residues","Relaxed Residues"],patch_artist=True,showmeans=True)
#colors=['pink','lightgreen']
#axis.set_ylim(-0.2,1.2)
#axis.set_xticklabels(["Strained Residues","Relaxed Residues"],fontsize=22)
#axis.tick_params(axis='y',which='major', labelsize=16)
#for element in ['boxes','whiskers', 'fliers', 'means', 'medians', 'caps']:
#    plt.setp(bplot[element], color='black')  
#for patch,color in zip(bplot['boxes'],colors):
#    patch.set_facecolor(color)
#axis.yaxis.grid(True)
#axis.set_ylabel("Normalized Fluctuations", fontsize=22)
    

'''
query="select bfactor from protein_strain_bfactors where molprob = 'Favored';"
cursor.execute(query)
results=cursor.fetchall()
unstrain_bfac_all=[]
for i in range(0,len(results)):
    unstrain_bfac_all.append(results[i][0])
    



stats.median_test(outlier_bfac,relax_bfactor,ties="above")
(385.310406004398, 8.686067719758959e-86, -0.0145, array([[  6310, 178494],
       [  8663, 176122]]))
stats.kstest(relax_bfactor,'norm')
KstestResult(statistic=0.1405527306444109, pvalue=0.0)
>>> stats.kstest(outlier_bfac,'norm')
KstestResult(statistic=0.11142764184625359, pvalue=6.678020884846427e-162)
stats.skew(relax_bfactor)
1.8120649040040155
>>> stats.skew(outlier_bfac)
2.1509526975275803'''

'''Outlier identification'''
#strain_q25=np.percentile(outlier_bfac,25)
#strain_q75=np.percentile(outlier_bfac,75)
#relax_q25=np.percentile(relax_bfactor,25)
#relax_q75=np.percentile(relax_bfactor,75)
#strain_iqr=strain_q75-strain_q25
#relax_iqr=relax_q75-relax_q25
#strain_cut_off=strain_iqr*3
#relax_cutoff=relax_iqr*3
#lower_strain_cutoff,upper_strain_cutoff=strain_q25-strain_cut_off,strain_q75+strain_cut_off
#lower_relax_cutoff,upper_relax_cutoff=relax_q25-relax_cutoff,relax_q75+relax_cutoff
#strain_filtered=[x for x in outlier_bfac if x>lower_strain_cutoff and x < upper_strain_cutoff]
#relax_filtered=[x for x in relax_bfactor if x>lower_relax_cutoff and x < upper_relax_cutoff]
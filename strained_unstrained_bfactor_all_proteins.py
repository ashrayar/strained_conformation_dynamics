# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 17:11:33 2019

@author: ashraya
"""

import MySQLdb
import numpy as np
from scipy import stats

#infile=open("/home/ashraya/MD_Project/outlier_molprob_analysis/strained_b-factor_all_proteins.csv","r")
#strained_bfactors=[]
#unstrained_bfactors=[]
#for line in infile:
#    lineparts=line.split(",")
#    bfactor=float(lineparts[4])
#    if bfactor<-9998:
#        continue
#    else:
#        strained_bfactors.append(bfactor)
#infile.close()
#
#db = MySQLdb.connect("localhost","root","admin","KINASE_PROJECT" )
#cursor = db.cursor() 
#infile=open("/home/ashraya/MD_Project/outlier_molprob_analysis/all_residues_bfactor.csv","r")
#for line in infile:
#    lineparts=line.split(",")
#    query="select bfactor from protein_strain_bfactors where pdb_id='"+lineparts[0]+"' and chain_id='"+lineparts[1]+"' and resno="+lineparts[2]+";"
#    x=cursor.execute(query)
#    results=cursor.fetchall()
#    if len(results)==0:
#        unstrained_bfactors.append(float(lineparts[-1][0:-1]))
#db.close()
#infile.close()

db = MySQLdb.connect("localhost","root","admin","KINASE_PROJECT" )
cursor = db.cursor() 
query="select bfactor from protein_strain_bfactors where molprob in ('Allowed','OUTLIER');"
x=cursor.execute(query)
strain_results=cursor.fetchall()
query="select bfactor from protein_strain_bfactors where molprob in ('Favored');"
x=cursor.execute(query)
unstrain_results=cursor.fetchall()
stat_results=stats.mannwhitneyu(strain_results,unstrain_results)
print np.mean(strain_results)
print np.mean(unstrain_results)
print stat_results
        
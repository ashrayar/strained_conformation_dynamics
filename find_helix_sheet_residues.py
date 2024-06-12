# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 23:57:05 2020

@author: ashraya
"""
'''
import MySQLdb
db = MySQLdb.connect("localhost","root","admin","KINASE_PROJECT" )
cursor = db.cursor() 

infile=open("/home/ashraya/MD_Project/all_protein_nma/max_flucts_raw_chains_termini_filtered.csv","r")
opfile1=open("/home/ashraya/MD_Project/all_protein_nma/helix_residues_bfactor.csv","w")
opfile2=open("/home/ashraya/MD_Project/all_protein_nma/sheet_residues_bfactor.csv","w")
for line1 in infile:
    if "nan\n" in line1:
        continue
    dsspfile=open("/home/ashraya/MD_Project/dssp_outputs/"+line1[0:4]+".dssp","r")
    helix=[]
    beta=[]
    flag=0
    for line in dsspfile:
        if line.startswith("  #  RESIDUE"):
            flag=1
            continue
        if flag==1:
            if "!*" in line:
                continue
            dssp=line[16]
            if dssp in ("E","B"):
                resno=line[5:10].strip()
                chain=line[11]
                beta.append(resno+","+chain)
            elif dssp == "H":
                resno=line[5:10].strip()
                chain=line[11]
                helix.append(resno+","+chain)
    dsspfile.close()
    query="select chain_id,resno,resname,bfactor,molprob from protein_strain_bfactors where pdb_id='"+line1[0:4]+"';"    
    cursor.execute(query)
    results=cursor.fetchall()
    for item in results:
        resno=str(item[1])
        to_search=resno+","+item[0]
        if to_search in helix:
            opfile1.write(line1[0:4]+","+item[0]+","+resno+","+item[2]+","+str(item[3])+","+item[4]+"\n")
        elif to_search in beta:
            opfile2.write(line1[0:4]+","+item[0]+","+resno+","+item[2]+","+str(item[3])+","+item[4]+"\n")
opfile1.close()
opfile2.close()
infile.close()'''

infile=open("/home/ashraya/MD_Project/outlier_molprob_analysis/molprob_outputs_list.csv","r")
relax_count=0
strain_count=0
for line1 in infile:
    if line1[0]=="#":
        continue
    print line1[0:-1]
    molprobfile=open("/home/ashraya/MD_Project/outlier_molprob_analysis/molprobity_outputs/"+line1[0:-1]+".rama","r")
    x=molprobfile.readline()
    strain_dict=[]
    relax_dict=[]
    for line in molprobfile:
        if line.startswith("SUMMARY"):
            continue
        lineparts=line.split(":")
        res_details=lineparts[0]
        chain=res_details[1]
        resno=res_details[2:6].strip()
        resname=res_details[6:].strip()
        if len(resname)>3:
            if resname[0]!="A":
                continue
            resname=resname[1:]
        if "Favored" in line:
            relax_dict.append(resno+"_"+chain)
        else:
            strain_dict.append(resno+"_"+chain)
    molprobfile.close()    
    dsspfile=open("/home/ashraya/MD_Project/dssp_outputs/"+line1[0:-1]+".dssp","r")
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
                if resno+"_"+chain in strain_dict:
                    strain_count+=1
                elif resno+"_"+chain in relax_dict:
                    relax_count+=1
    dsspfile.close()
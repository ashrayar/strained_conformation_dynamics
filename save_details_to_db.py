# -*- coding: utf-8 -*-
"""
Created on Thu May 30 16:01:01 2019

@author: ashraya
"""
'''
/home/ashraya/MD_Project/hrd_molprob_all.csv
/home/ashraya/MD_Project/kinase_list_latest_with_resolution.csv
/home/ashraya/MD_Project/pdb_uniprot_name_mapping_all.csv
/home/ashraya/MD_Project/all_kinase_all_locations.csv
/home/ashraya/MD_Project/dfg_molprob_all.csv
/home/ashraya/MD_Project/all_kinase_active_final.csv
/home/ashraya/MD_Project/all_kinase_inactive.csv
/home/ashraya/MD_Project/hrd_flucts_strained_all.csv
/home/ashraya/MD_Project/hrd_flucts_unstrained_all.csv
/home/ashraya/MD_Project/dfg_flucts_unstrained_all.csv
/home/ashraya/MD_Project/dfg_flucts_strained_all.csv
/home/ashraya/MD_Project/strained_dfg_b-factor_all.csv
/home/ashraya/MD_Project/strained_hrd_b-factor_all.csv
/home/ashraya/MD_Project/unstrained_dfg_b-factor_all.csv
/home/ashraya/MD_Project/unstrained_hrd_b-factor_all.csv
'''

kinases=dict()
list_of_pdbs=[]
infile=open("/home/ashraya/MD_Project/lower_res_analysis/1.7_to_2_all_kinases_resolution.csv","r")
for line in infile:
    lineparts=line.split(",")
    pdb_id=lineparts[0]
    list_of_pdbs.append(pdb_id)
infile.close()

'''resolution'''
for item in list_of_pdbs:
    #pdb_id=item[0:4].upper()
    infile=open("lower_res_analysis/1.7_to_2_all_kinases_resolution.csv","r")
    #x=infile.readline()
    for line in infile:
        if line[0:6]==pdb_id:
            lineparts=line.split(",")
            resol=float(lineparts[1][0:-1])
            kinases[item]=[resol]
            break
    infile.close()

'''uniprot mapping'''
infile=open("/home/ashraya/MD_Project/lower_res_analysis/1.7_to_2_ser_thr_kinases_pfam_mapping.csv","r")
for line in infile:
    lineparts=line.split(",")
    pdb_id=lineparts[3].lower()+"_"+lineparts[4]
    kinases[pdb_id].append(lineparts[0])
infile.close()

infile=open("/home/ashraya/MD_Project/lower_res_analysis/1.7_to_2_tyr_kinases_pfam_mapping.csv","r")
for line in infile:
    lineparts=line.split(",")
    pdb_id=lineparts[3].lower()+"_"+lineparts[4]
    kinases[pdb_id].append(lineparts[0])
infile.close()

'''hrd arg resno, phi, psi, strain'''
infile=open("/home/ashraya/MD_Project/lower_res_analysis/all_kinases_hrd_molprob.csv","r")
for line in infile:
    lineparts=line.split(",")
    arg_resno=int(lineparts[1])
    arg_phi=int(lineparts[3])
    arg_psi=int(lineparts[4])
    if lineparts[2] == "AARG":
        resname="ARG"
    else:
        resname=lineparts[2]
    strain=lineparts[-1][0:-1]
    kinases[lineparts[0]].extend((arg_resno,resname,arg_phi,arg_psi,strain))
infile.close()

'''other residue locations'''
infile=open("/home/ashraya/MD_Project/lower_res_analysis/low_res_all_locations.csv","r")
x=infile.readline()
for line in infile:
    lineparts=line.split(",")
    dfg_asp=int(lineparts[4])
    ehelix_his=int(lineparts[2])
    ehelix_his_residue=lineparts[5]
    fhelix_asp=int(lineparts[3])
    glu=int(lineparts[7][0:-1])
    lys=int(lineparts[6])
    kinases[lineparts[0]].extend((dfg_asp,ehelix_his,ehelix_his_residue,fhelix_asp,lys,glu))
infile.close()

'''dfg asp phi, psi, strain'''
'''for item in list_of_pdbs:
    flag=0
    infile=open("dfg_molprob_all.csv","r")
    for line in infile:
        if line[0:6]==item:
            lineparts=line.split(",")
            dfg_phi=int(lineparts[3])
            dfg_psi=int(lineparts[4])
            dfg_strain=lineparts[-1][0:-1]
            flag=1
            kinases[lineparts[0]].extend((dfg_phi,dfg_psi,dfg_strain))
            break
    infile.close()
    if flag==0:
        kinases[item].extend((-999,-999,"missing"))'''
        
infile=open("/home/ashraya/MD_Project/lower_res_analysis/all_kinases_dfg_molprob.csv","r")
for line in infile:
    lineparts=line.split(",")
    arg_resno=int(lineparts[1])
    arg_phi=int(lineparts[3])
    arg_psi=int(lineparts[4])
    if lineparts[2] == "AASP":
        resname="ASP"
    else:
        resname=lineparts[2]
    strain=lineparts[-1][0:-1]
    kinases[lineparts[0]].extend((arg_phi,arg_psi,strain))
infile.close()

'''kinase activity'''
#active=[]
#infile=open("all_kinase_active_final.csv","r")
#for line in infile:
#    active.append(line[0:6])
#infile.close()
#for item in list_of_pdbs:
#    if item in active:
#        kinases[item].append("active")
#    else:
#        kinases[item].append("inactive")

for item in list_of_pdbs:
    kinases[item].append('unknown')

'''hrd flucts'''
for item in list_of_pdbs:
    flag=0
    #infile=open("/home/ashraya/MD_Project/hrd_flucts_all.csv","r")
    infile=open("/home/ashraya/MD_Project/lower_res_analysis/all_kinases_hrd_flucts.csv","r")
    for line in infile:
        if line[0:6]==item:
            lineparts=line.split(",")
            hrd_fluct=float(lineparts[-1][0:-1])
            flag=1
            kinases[lineparts[0]].append(hrd_fluct)
            break
    infile.close()
    if flag==0:
        kinases[item].append(-999.0)


'''dfg flucts'''
for item in list_of_pdbs:
    flag=0
    infile=open("/home/ashraya/MD_Project/lower_res_analysis/all_kinases_dfg_flucts.csv","r")
    for line in infile:
        if line[0:6]==item:
            lineparts=line.split(",")
            dfg_fluct=float(lineparts[-1][0:-1])
            flag=1
            kinases[lineparts[0]].append(dfg_fluct)
            break
    infile.close()
    if flag==0:
        kinases[item].append(-999.0)
            
'''hrd bfactor'''
#infile=open("/home/ashraya/MD_Project/strained_hrd_b-factor_all.csv","r")
#for line in infile:
#    lineparts=line.split(",")
#    bfac=float(lineparts[-1][0:-1])
#    kinases[lineparts[0]].append(bfac)
#infile.close()
#infile=open("/home/ashraya/MD_Project/unstrained_hrd_b-factor_all.csv","r")
#for line in infile:
#    lineparts=line.split(",")
#    bfac=float(lineparts[-1][0:-1])
#    kinases[lineparts[0]].append(bfac)
#infile.close()

infile=open("/home/ashraya/MD_Project/lower_res_analysis/hrd_b-factor_all.csv","r")
for line in infile:
    lineparts=line.split(",")
    bfac=float(lineparts[-1][0:-1])
    kinases[lineparts[0]].append(bfac)
infile.close()

'''dfg bfactor'''
for item in list_of_pdbs:
    flag=0
    infile=open("/home/ashraya/MD_Project/lower_res_analysis/dfg_b-factor_all.csv","r")
    for line in infile:
        if line[0:6]==item:
            lineparts=line.split(",")
            dfg_fluct=float(lineparts[-1][0:-1])
            flag=1
            kinases[lineparts[0]].append(dfg_fluct)
            break
    infile.close()
    if flag==0:
        kinases[item].append(-999.0)

infile=open("/home/ashraya/MD_Project/lower_res_analysis/pfam_domain_kinase_family_mapping.csv","r")
uni_fam_map=dict()
for line in infile:
    lineparts=line.split(",")
    uni_fam_map[lineparts[0]]=lineparts[1][0:-1]
infile.close()
for item in list_of_pdbs:
    kinases[item].append(uni_fam_map[kinases[item][1]])

import MySQLdb
db = MySQLdb.connect("localhost","root","admin","KINASE_PROJECT" )
cursor = db.cursor() 
for item in list_of_pdbs:
    insertstr="Insert into kinases values('%s','%f','%s','%d','%s','%d','%d','%s','%d','%d','%d','%d','%d','%d','%d','%s','%s','%s','%f','%f','%f','%f','%f','%f','%s')" % \
    (item,kinases[item][0],kinases[item][1],kinases[item][2],kinases[item][3],kinases[item][7],kinases[item][8],kinases[item][9],kinases[item][10],kinases[item][11],kinases[item][12],kinases[item][4],kinases[item][5],kinases[item][13],kinases[item][14],kinases[item][16],kinases[item][6],kinases[item][15],-999.0,kinases[item][17],-999.0,kinases[item][18],kinases[item][19],kinases[item][20],kinases[item][21])
    #(item,kinases[item][0],kinases[item][1],kinases[item][2],kinases[item][3],kinases[item][7],kinases[item][8],kinases[item][9],kinases[item][10],kinases[item][11],kinases[item][12],kinases[item][4],kinases[item][5],kinases[item][13],kinases[item][14],kinases[item][16],kinases[item][6],kinases[item][15],kinases[item][17],kinases[item][18],kinases[item][19],kinases[item][20],kinases[item][21])
    
    cursor.execute(insertstr)
    db.commit()
db.close()


db = MySQLdb.connect("localhost","root","admin","KINASE_PROJECT" )
cursor = db.cursor() 
for item in list_of_pdbs:
    infile=open("/home/ashraya/MD_Project/lower_res_analysis/h_bond_counts_all.csv","r")
    x=infile.readline()
    flag=0
    for line in infile:
        lineparts=line.split(",")
        if lineparts[0]==item:
            insertstr="insert into kinase_hbonds values('%s','%d','%d','%d','%d','%d','%d')" % \
            (item,int(lineparts[1]),int(lineparts[2]),int(lineparts[3]),int(lineparts[4]),int(lineparts[5]),int(lineparts[6]))
            cursor.execute(insertstr)
            db.commit()
            flag=1
            break
    infile.close()
    if flag==0:
        insertstr="insert into kinase_hbonds values('%s','%d','%d','%d','%d','%d','%d')" % \
        (item,0,0,0,0,0,0)
        cursor.execute(insertstr)
        db.commit()
db.close()


for item in list_of_pdbs:
    flag=0            
    infile=open("/home/ashraya/MD_Project/dfg_flucts_no_ligands_all.csv","r")
    for line in infile:
        if line[0:6]==item:
            lineparts=line.split(",")
            hrd_fluct=float(lineparts[-1][0:-1])
            flag=1
            insertstr="update kinases set dfg_fluct_no_ligand='%f' where pdb_id='%s';" %\
            (hrd_fluct,lineparts[0])
            cursor.execute(insertstr)
            db.commit()
            break
    infile.close()
    if flag==0:
        insertstr="update kinases set dfg_fluct_no_ligand='%f' where pdb_id='%s';" %\
        (-999.0,item)
        cursor.execute(insertstr)
        db.commit()
db.close()

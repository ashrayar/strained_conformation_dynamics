# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 18:17:24 2019

@author: ashraya
"""
'''
nma_structs=dict()
exclude=["1k2a_A","1mwq_A","1mwq_B","3ll1_A","3ll2_A","5jsk_B"]
infile=open("/home/ashraya/MD_Project/all_protein_nma/structures_for_nma.csv","r")
for line in infile:
    lineparts=line.split(",")
    if lineparts[0] not in exclude:
        nma_structs[lineparts[0]]=[int(lineparts[1]),int(lineparts[2][0:-1])]
infile.close()

structnames=nma_structs.keys()

infile=open("/home/ashraya/MD_Project/loop_residue_bfactors_1A.csv","r")
opfile=open("/home/ashraya/MD_Project/all_protein_nma/nma_loop_flucts.csv","w")
prev_struct=" "
loop_res=dict()
for line in infile:
    lineparts=line.split(",")
    struct=lineparts[0]+"_"+lineparts[1]
    if prev_struct!=struct:
        if prev_struct in structnames:
            fluctfile=open("/home/ashraya/MD_Project/all_protein_nma/fluctuations/"+prev_struct+"_flucts.csv","r")
            for line1 in fluctfile:
                line1parts=line1.split("_")
                res=int(line1parts[0])
                fluctuation=line1parts[2].split(",")[1][0:-1]
                loop_res_res=loop_res.keys()
                if res in loop_res_res:
                    opfile.write(prev_struct+","+str(res)+","+loop_res[res][0]+","+fluctuation+","+loop_res[res][1]+"\n")
            fluctfile.close()
        loop_res=dict()
    if struct in structnames:
        resno=int(lineparts[2])
        fluct_res=resno-nma_structs[struct][0]+1
        loop_res[fluct_res]=[lineparts[3],lineparts[5][0:-1]]
    prev_struct=struct
infile.close()
opfile.close()'''

nma_structs=[]
#infile=open("/home/ashraya/MD_Project/no_missing_all_proteins_new.csv","r")
infile=open("/home/ashraya/MD_Project/all_protein_nma/max_flucts_raw_chains_termini_filtered.csv","r")
for line in infile:
    #nma_structs.append(line[0:-1])
    if "nan\n" in line:
        continue
    nma_structs.append(line[0:4])
infile.close()
#exclude=["1knl","1m1n"]
exclude=[]
infile=open("/home/ashraya/MD_Project/loop_residue_bfactors.csv","r")
#opfile=open("/home/ashraya/MD_Project/all_protein_nma/nma_loop_flucts_1.2.csv","w")
opfile=open("/home/ashraya/MD_Project/all_protein_nma/nma_loop_flucts_all_norm.csv","w")
for line in infile:
    if line[0:4] in exclude:
        continue
    if line[0:4] in nma_structs:
        lineparts=line[0:-1].split(",")
        res=lineparts[2]
        chain=lineparts[1]
        #fluctfile=open("/home/ashraya/MD_Project/all_protein_nma/fluctuations_new/"+line[0:4]+"_flucts.csv","r")
        fluctfile=open("/home/ashraya/MD_Project/all_protein_nma/fluctuations_raw/normalized_flucts/"+line[0:4]+"_norm_flucts.csv","r")
        for line1 in fluctfile:
            line1parts=line1[0:-1].split("_")
            if line1parts[1]==chain and line1parts[0]==res:
                #opfile.write(lineparts[0]+","+lineparts[1]+","+res+","+lineparts[3]+","+line1[0:-1].split(",")[1]+","+lineparts[5]+"\n")
                opfile.write(lineparts[0]+","+lineparts[1]+","+res+","+lineparts[3]+","+line1[0:-1].split(",")[1]+","+line1[0:-1].split(",")[2]+","+lineparts[5]+"\n")
                break
        fluctfile.close()
infile.close()
opfile.close()
                
            
            
            


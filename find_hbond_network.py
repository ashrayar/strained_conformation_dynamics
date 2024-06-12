# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 11:54:35 2019

@author: ashraya
"""

infile=open("/home/ashraya/MD_Project/lower_res_analysis/low_res_all_locations.csv","r")
opfile=open("/home/ashraya/MD_Project/lower_res_analysis/h_bonds_detected_all.csv","w")
x=infile.readline()
for line1 in infile:
    line1parts=line1.split(",")
    ehis=int(line1parts[2])
    hrd_his=int(line1parts[1])
    f_asp=int(line1parts[3])
    dfg_asp=int(line1parts[4])
    pdb_id=line1parts[0]
    if ehis==-999 or hrd_his==-999 or f_asp==-999 or dfg_asp==-999:
        continue
    else:
        ion_flag=0
        ms_flag=0
        int_file=open("/home/ashraya/MD_Project/lower_res_analysis/pic_input_and_output/"+pdb_id+".int","r")        
        for line in int_file:
            if line.startswith(" Intraprotein Ionic Interaction"):
#                print "ion flag set 1"
                ion_flag=1
                continue
            if ion_flag==1:
                ion_flag=2
                continue
            if ion_flag==2:
                lineparts=line.split("\t")
                if len(lineparts)>1:
                    if int(lineparts[0])==ehis:
                        if int(lineparts[3])==f_asp:
                            opfile.write(pdb_id+",ionic,E-His,F-Asp,"+str(ehis)+","+str(f_asp)+"\n")
                            ion_flag=0
                            continue
                if line.startswith(" Intraprotein aromatic-aromatic Interaction"):
                    ion_flag=0
            if line.startswith(" Intraprotein Main chain-Side chain"):
                ms_flag=1
                print "ms flag set 1"
                continue
            if line.startswith("POS\tCHAIN") and ms_flag==1:
                ms_flag=2
                print "ms flag set 2"
                continue
            if ms_flag==2:
                ms_flag=3
                continue
            if ms_flag==3:
#                print "ms flag 3"
                lineparts=line.split("\t")
                if lineparts[0]=="\n":
                    continue
                if line.startswith(" Intraprotein Side chain-Side chain"):
                    ms_flag=0
                    break
                if int(lineparts[0])==hrd_his:
                    if int(lineparts[4])==dfg_asp:
                        if lineparts[3].strip()=="NE2" and lineparts[7].strip()=="O":
                            opfile.write(pdb_id+",MC-SC,HRD-His,DFG-Asp,"+str(hrd_his)+","+lineparts[3].strip()+","+str(dfg_asp)+","+lineparts[7].strip()+"\n")
                    elif int(lineparts[4])==f_asp:
                        if lineparts[3].strip()=="N" and (lineparts[7].strip()=="OD1" or lineparts[7].strip()=="OD2"):
                            opfile.write(pdb_id+",MC-SC,HRD-His,F-Asp,"+str(hrd_his)+","+lineparts[3].strip()+","+str(f_asp)+","+lineparts[7].strip()+"\n")
                elif int(lineparts[0])==hrd_his+1:
                    if int(lineparts[4])==f_asp:
                        if lineparts[3].strip()=="N" and (lineparts[7].strip()=="OD1" or lineparts[7].strip()=="OD2"):
                            opfile.write(pdb_id+",MC-SC,HRD-Arg,F-Asp,"+str(hrd_his+1)+","+lineparts[3].strip()+","+str(f_asp)+","+lineparts[7].strip()+"\n")
                elif int(lineparts[0])==hrd_his+2:
                    if int(lineparts[4])==hrd_his:
                        if lineparts[3].strip()=="N" and lineparts[7].strip()=="ND1":
                            opfile.write(pdb_id+",MC-SC,HRD-Asp,HRD-His,"+str(hrd_his+2)+","+lineparts[3].strip()+","+str(hrd_his)+","+lineparts[7].strip()+"\n")
        int_file.close()
infile.close()
opfile.close()
            
helplibrary(bio3d)
kin_strained<-read.csv("/home/ashraya/MD_Project/strained_no_missing.txt",header = FALSE)
kin_strained<-kin_strained$V1
for(j in 1:length(kin_strained))
{
  filename=paste("/home/ashraya/Downloads/corrected_structs/",kin_strained[j],"_fixed.pdb",sep="")
  kin_pdb<-read.pdb(filename)
  aanm<-aanma.pdb(kin_pdb,keep=17,outmodes = "calpha")
  flucts_norm<-aanm$fluctuations/max(aanm$fluctuations)
  file1name=paste("/home/ashraya/MD_Project/fluctuations/",kin_strained[j],"_strained_flucts.csv",sep="")
  for(k in 1:length(flucts_norm))
  {
    writestr=paste(names(aanm$fluctuations[k]),flucts_norm[k],sep = ",")
    cat(writestr,file = file1name,append = TRUE,sep = "\n")
  }
}
kin_unstrained<-read.csv("/home/ashraya/MD_Project/unstrained_no_missing.txt",header = FALSE)
kin_unstrained<-kin_unstrained$V1
for(j in 22:length(kin_unstrained))
{
  print(kin_unstrained[j])
  filename=paste("/home/ashraya/Downloads/corrected_structs/",kin_unstrained[j],"_fixed.pdb",sep="")
  kin_pdb<-read.pdb(filename)
  aanm<-aanma.pdb(kin_pdb,keep=17,outmodes = "calpha")
  flucts_norm<-aanm$fluctuations/max(aanm$fluctuations)
  file1name=paste("/home/ashraya/MD_Project/fluctuations/",kin_unstrained[j],"_unstrained_flucts.csv",sep="")
  for(k in 1:length(flucts_norm))
  {
    writestr=paste(names(aanm$fluctuations[k]),flucts_norm[k],sep = ",")
    cat(writestr,file = file1name,append = TRUE,sep = "\n")
  }
}

kin_strained<-read.csv("/home/ashraya/MD_Project/missing_residue_details_temp_remaining.txt",header = FALSE)
kin_strained<-kin_strained$V1
kin_strained<-c("5n3r_A")
for(j in 1:length(kin_strained))
{
  filename=paste("/home/ashraya/MD_Project/corrected_structs/",kin_strained[j],"_fixed.pdb",sep="")
  kin_pdb<-read.pdb(filename)
  aanm<-aanma.pdb(kin_pdb,keep=17,outmodes = "calpha")
  flucts_norm<-aanm$fluctuations/max(aanm$fluctuations)
  file1name=paste("/home/ashraya/MD_Project/fluctuations_new/",kin_strained[j],"_flucts.csv",sep="")
  for(k in 1:length(flucts_norm))
  {
    writestr=paste(names(aanm$fluctuations[k]),flucts_norm[k],sep = ",")
    cat(writestr,file = file1name,append = TRUE,sep = "\n")
  }
}
kin_unstrained<-c("5n1e_A")
kin_unstrained<-read.csv("/home/ashraya/MD_Project/unstrained_no_missing_all.csv",header = FALSE)
kin_unstrained<-kin_unstrained$V1
for(j in 1:length(kin_unstrained))
{
  print(kin_unstrained[j])
  filename=paste("/home/ashraya/MD_Project/corrected_structs/",kin_unstrained[j],"_fixed.pdb",sep="")
  kin_pdb<-read.pdb(filename)
  aanm<-aanma.pdb(kin_pdb,keep=17,outmodes = "calpha",ligand=FALSE)
  flucts_norm<-aanm$fluctuations/max(aanm$fluctuations)
  file1name=paste("/home/ashraya/MD_Project/",kin_unstrained[j],"_unstrained_flucts_test_2.csv",sep="")
  for(k in 1:length(flucts_norm))
  {
    writestr=paste(names(aanm$fluctuations[k]),flucts_norm[k],sep = ",")
    cat(writestr,file = file1name,append = TRUE,sep = "\n")
  }
}
filenames=c("/home/ashraya/MD_Project/kinase_structures_renumbered/1gz8_A_renumbered.pdb","/home/ashraya/MD_Project/kinase_structures_renumbered/1fmk_A_renumbered.pdb","/home/ashraya/MD_Project/kinase_structures_renumbered/1jks_A_renumbered.pdb","/home/ashraya/MD_Project/kinase_structures_renumbered/1p4o_A_renumbered.pdb","/home/ashraya/MD_Project/kinase_structures_renumbered/1p4o_B_renumbered.pdb","/home/ashraya/MD_Project/kinase_structures_renumbered/1rdq_E_renumbered.pdb")

kin_unstrained<-read.csv("/home/ashraya/MD_Project/hrd_flucts_all.csv",header = FALSE)
kin_unstrained<-kin_unstrained$V1
for(j in 1:length(kin_unstrained))
{
  print(kin_unstrained[j])
  filename=paste("/home/ashraya/MD_Project/corrected_structs_new/",kin_unstrained[j],"_fixed_new.pdb",sep="")
  kin_pdb<-read.pdb(filename)
  aanm<-aanma.pdb(kin_pdb,keep=17,outmodes = "calpha")
  flucts_norm<-aanm$fluctuations/max(aanm$fluctuations)
  file1name=paste("/home/ashraya/MD_Project/fluctuations_no_ligands/",kin_unstrained[j],"_flucts_no_ligands.csv",sep="")
  for(k in 1:length(flucts_norm))
  {
    writestr=paste(names(aanm$fluctuations[k]),flucts_norm[k],sep = ",")
    cat(writestr,file = file1name,append = TRUE,sep = "\n")
  }
}


kin_all<-read.csv("/home/ashraya/MD_Project/lower_res_analysis/nma_pdb_remaining.csv",header = FALSE)
kin_all<-kin_all$V1
for(j in 1:length(kin_all))
{
  print(kin_all[j])
  filename=paste("/home/ashraya/MD_Project/lower_res_analysis/corrected_structs/",kin_all[j],"_fixed.pdb",sep="")
  kin_pdb<-read.pdb(filename)
  aanm<-aanma.pdb(kin_pdb,keep=17,outmodes = "calpha")
  flucts_norm<-aanm$fluctuations/max(aanm$fluctuations)
  file1name=paste("/home/ashraya/MD_Project/lower_res_analysis/fluctuations/",kin_all[j],"_flucts.csv",sep="")
  for(k in 1:length(flucts_norm))
  {
    writestr=paste(names(aanm$fluctuations[k]),flucts_norm[k],sep = ",")
    cat(writestr,file = file1name,append = TRUE,sep = "\n")
  }
}

all_proteins<-read.csv("/home/ashraya/MD_Project/all_protein_nma/structures_for_nma.csv",header=FALSE)
all_proteins<-all_proteins$V1
for(j in 1:length(all_proteins))
{
  print(all_proteins[j])
  filename=paste("/home/ashraya/MD_Project/all_protein_nma/corrected_structs/",all_proteins[j],"_fixed.pdb",sep="")
  kin_pdb<-read.pdb(filename)
  aanm<-aanma.pdb(kin_pdb,keep=17,outmodes = "calpha")
  flucts_norm<-aanm$fluctuations/max(aanm$fluctuations)
  file1name=paste("/home/ashraya/MD_Project/all_protein_nma/fluctuations/",all_proteins[j],"_flucts.csv",sep="")
  for(k in 1:length(flucts_norm))
  {
    writestr=paste(names(aanm$fluctuations[k]),flucts_norm[k],sep = ",")
    cat(writestr,file = file1name,append = TRUE,sep = "\n")
  }
}


all_proteins<-read.csv("/home/ashraya/MD_Project/all_protein_nma/flucts_raw_to_calculate.csv",header=FALSE)
all_proteins<-all_proteins$V1
for(j in 18:length(all_proteins))
{
  print(all_proteins[j])
  filename=paste("/home/ashraya/MD_Project/all_protein_nma/corrected_structs_new/",all_proteins[j],"_fixed.pdb",sep="")
  kin_pdb<-read.pdb(filename)
  aanm<-aanma.pdb(kin_pdb,keep=17,outmodes = "calpha")
  #flucts_norm<-aanm$fluctuations/max(aanm$fluctuations)
  #file1name=paste("/home/ashraya/MD_Project/all_protein_nma/fluctuations_new/",all_proteins[j],"_flucts.csv",sep="")
  file1name=paste("/home/ashraya/MD_Project/all_protein_nma/fluctuations_raw/",all_proteins[j],"_flucts.csv",sep="")
  #for(k in 1:length(flucts_norm))
  #{
    #writestr=paste(names(aanm$fluctuations[k]),flucts_norm[k],sep = ",")
    #cat(writestr,file = file1name,append = TRUE,sep = "\n")
  #}
  for(k in 1:length(aanm$fluctuations))
  {
    writestr=paste(names(aanm$fluctuations[k]),aanm$fluctuations[k],sep = ",")
    cat(writestr,file = file1name,append = TRUE,sep = "\n")
  }
}

kin_all<-read.csv("/home/ashraya/MD_Project/kinase_nma_list.csv",header=FALSE)
kin_all<-kin_all$V1
for(j in 1:length(kin_all))
{
  print(kin_all[j])
  filename=paste("/home/ashraya/MD_Project/corrected_structs_new/",kin_all[j],"_fixed_new.pdb",sep="")
  kin_pdb<-read.pdb(filename)
  aanm<-aanma.pdb(kin_pdb,keep=17,outmodes = "calpha")
  file1name=paste("/home/ashraya/MD_Project/fluctuations_no_ligands/fluctuations_raw/",kin_all[j],"_flucts.csv",sep="")
  for(k in 1:length(aanm$fluctuations))
  {
    writestr=paste(names(aanm$fluctuations[k]),aanm$fluctuations[k],sep = ",")
    cat(writestr,file = file1name,append = TRUE,sep = "\n")
  }
}
pdb_1om1<-read.pdb("/home/ashraya/MD_Project/kinase_structures_renumbered/1gz8_A_renumbered.pdb")
pdb_1om1<-trim.pdb(pdb_1om1,atom.select(pdb_1om1,"protein",chain="A"))
aanm_1om1<-aanma.pdb(pdb_1om1,keep=7,outmodes = "calpha")
flucts_1om1_norm<-aanm_1om1$fluctuations/max(aanm_1om1$fluctuations)
filecon<-file("1gz8_flucts_new.csv")
for(j in 1:length(flucts_1om1_norm))
{
  writestr=paste(names(aanm_1om1$fluctuations[j]),flucts_1om1_norm[j],sep = ",")
  cat(writestr,file = "1gz8_flucts_new.csv",append = TRUE,sep = "\n")
  
}
close(filecon)

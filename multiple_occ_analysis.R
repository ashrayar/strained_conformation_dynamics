library(bio3d)

list_of_pdbs<-read.csv("/home/ashraya/MD_Project/multiple_occ_analysis/structures_1.0_june_2019.txt",header=FALSE)
pdblist<-list_of_pdbs$V1
pdblist<-c("2ewk")
for(j in 103:length(pdblist))
{
  print(j)
  filename=paste("/home/ashraya/Refinement_project/trial_2/structures_1.2/",pdblist[j],".pdb",sep="")
  test_pdb<-read.pdb(filename,rm.alt = FALSE)
  alt_array=test_pdb$atom$alt
  a_index<-c()
  b_index<-c()
  for(i in 1:length(test_pdb$atom$alt))
  {
    if (!(is.na(test_pdb$atom$alt[i])))
    {
      if(test_pdb$atom$alt[i]=="A")
      {
        a_index<-c(a_index,test_pdb$atom$eleno[i])
      }
      if(test_pdb$atom$alt[i]=="B")
      {
        b_index<-c(b_index,test_pdb$atom$eleno[i])
      }
    }
  }
  if(length(a_index)==0)
  {
    a_index=which(alt_array=="A")
  }
  if(length(b_index)==0)
  {
    b_index=which(alt_array=="B")
  }
  cut_pdb_a<-atom.select(test_pdb,eleno = b_index,inverse = TRUE)
  trunc_pdb_a<-trim.pdb(test_pdb,cut_pdb_a)
  writefile_a=paste("/home/ashraya/MD_Project/multiple_occ_analysis/split_pdb/",pdblist[j],"_A_occ.pdb",sep="")
  write.pdb(trunc_pdb_a,file=writefile_a)
  cut_pdb_b<-atom.select(test_pdb,eleno = a_index,inverse = TRUE)
  trunc_pdb_b<-trim.pdb(test_pdb,cut_pdb_b)
  writefile_b=paste("/home/ashraya/MD_Project/multiple_occ_analysis/split_pdb/",pdblist[j],"_B_occ.pdb",sep="")
  write.pdb(trunc_pdb_b,file=writefile_b)

}

filename=paste("/home/ashraya/Refinement_project/trial_2/structures_1.2/2ce2.pdb",sep="")
test_pdb<-read.pdb(filename,rm.alt = FALSE)
alt_array=test_pdb$atom$alt
a_index<-c()
b_index<-c()
for(i in 1:length(test_pdb$atom$alt))
{
  if (!(is.na(test_pdb$atom$alt[i])))
  {
    if(test_pdb$atom$alt[i]=="A")
    {
      a_index<-c(a_index,test_pdb$atom$eleno[i])
    }
    if(test_pdb$atom$alt[i]=="B")
    {
      b_index<-c(b_index,test_pdb$atom$eleno[i])
    }
  }
}
cut_pdb_a<-atom.select(test_pdb,eleno = b_index,inverse = TRUE)
trunc_pdb_a<-trim.pdb(test_pdb,cut_pdb_a)
writefile_a=paste("/home/ashraya/MD_Project/multiple_occ_analysis/split_pdb/2ce2_A_occ.pdb",sep="")
write.pdb(trunc_pdb_a,file=writefile_a)
cut_pdb_b<-atom.select(test_pdb,eleno = a_index,inverse = TRUE)
trunc_pdb_b<-trim.pdb(test_pdb,cut_pdb_b)
writefile_b=paste("/home/ashraya/MD_Project/multiple_occ_analysis/split_pdb/2ce2_B_occ.pdb",sep="")
write.pdb(trunc_pdb_b,file=writefile_b)

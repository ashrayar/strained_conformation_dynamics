library(bio3d)
kin_strained<-read.csv("/home/ashraya/MD_Project/kinase_strained.csv",header = FALSE)
kin_strained<-kin_strained$V1
strain_pdbs<-vector(mode="character",length=60)
for(j in 1:length(kin_strained))
{
  filename=paste("/home/ashraya/MD_Project/kinase_structures_renumbered/",kin_strained[j],"_renumbered.pdb",sep="")
  strain_pdbs[j]<-filename
}
aln_new <- pdbaln(strain_pdbs)
pdbs_new <- read.all(aln_new)
modes_strained<-aanma.pdbs(pdbs_new,subspace = 2)

#trial NMA
library(bio3d)
#list_file<-read.csv("/home/ashraya/MD_Project/pfam_mapping_latest_high_res_only.csv",header=FALSE)
list_file<-read.csv("/home/ashraya/MD_Project/lower_res_analysis/pdb_range_all_kinases_final.csv",header=FALSE)
colnames(list_file)=c("pdb_id","start","end")
pdb_list<-list_file$pdb_id
pdb_list<-list.files("/home/ashraya/MD_Project/lower_res_analysis/kinase_structures_truncated")
#pdb_list=c("3ckw_A")
for(j in 1:length(pdb_list)){
  pdbfile=paste("/home/ashraya/MD_Project/lower_res_analysis/truncated_structs_new/",pdb_list[j],"_trunc_new.pdb",sep="")
  pdb_trunc=read.pdb(pdbfile)
  #lineparts=strsplit(pdb_list[j],"_")[[1]]
  #writefile=paste("/home/ashraya/MD_Project/kinase_structures_renumbered/",lineparts[1],"_",lineparts[2],"_renumbered.pdb",sep="")
  writefile=paste("/home/ashraya/MD_Project/lower_res_analysis/kinase_structures_renumbered/",pdb_list[j],"_renumbered.pdb",sep="")
  write.pdb(pdb_trunc,file=writefile,resno = pdb_trunc$atom$resno-pdb_trunc$atom$resno[1]+1)
}

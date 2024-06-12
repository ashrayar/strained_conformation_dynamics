ids<-read.csv("/home/ashraya/MD_Project/outlier_molprob_analysis/cluster_12_sample_pca_structs.csv",header=FALSE)
ids<-ids$V1
#cluster 3
ids<-c("1LUG_A","2EU2_A","2FOS_A","2FOU_A","2FOV_A","2ILI_A","2WEG_A","3D92_A","3D93_A","3HS4_A","3KS3_A","3M2Y_A","3MHO_A","3RJ7_A","3SAX_A","3V7X_A","3VBD_A","4FPT_A","4FRC_A","4FU5_A","4FVN_A","4FVO_A","4ITO_A","4MTY_A","4PZH_A","4Q6D_A","4Q6E_A","4Q78_A","4RUX_A","4RUY_A","4WL4_A","4WW6_A","4YX4_A","4YXI_A","4YXO_A","4YXU_A","4YYT_A","5BYI_A","5DOH_A","5DOH_B","5DRS_A","5LJQ_A","5LJT_A","5LL4_A","5LL4_B","5LL8_A","5LLC_A","5LLG_A","5M78_A","5MJN_A","5NXG_A","5NXO_A","5NXV_A","5NXW_A","5NY1_A","5NY6_A","5NYA_A","5Y2R_A","5Y2S_A","5YUI_A","6B4D_A","6G3Q_A","6G6T_A","6GDC_A","6GM9_A","6IC2_A")
raw.files<-get.pdb(ids)
files <- pdbsplit(raw.files, ids)
pdbs <- pdbaln(files)
pdbs$id <- substr(basename(pdbs$id),1,6)
seqidentity(pdbs)
rmsd(pdbs,fit=TRUE)
xyz=pdbfit(pdbs)
rd <- rmsd(xyz)
hist(rd, breaks=40, xlab="RMSD (Å)", main="Histogram of RMSD")
hc.rd <- hclust(as.dist(rd))
pdbs$id <- substr(basename(pdbs$id), 1, 6)
hclustplot(hc.rd,k=3,labels=pdbs$id, cex=0.5,ylab="RMSD (Å)", main="RMSD Cluster Dendrogram", fillbox=FALSE)
pc.xray <- pca.xyz(xyz,rm.gaps = TRUE)
plot(pc.xray)
plot(pc.xray, pc.axes=1:2)
identify(pc.xray$z[,1:2], labels=basename.pdb(pdbs$id))
par(mfrow = c(2, 1), cex = 0.75, mar = c(3, 4, 1, 1))
plot.bio3d(pc.xray$au[,1], resno= pdbs$resno[2,],sse=pdbs$sse[2,],ylab="PC1")
plot.bio3d(pc.xray$au[,2], ylab="PC2")
file1name="/home/ashraya/MD_Project/cluster_3_pca_residue_contribution.csv"
for(k in 1:length(pc.xray$au[,1]))
{
  writestr=paste(k,pc.xray$au[k,1],pc.xray$au[k,2],sep = ",")
  cat(writestr,file = file1name,append = TRUE,sep = "\n")
}
file1name="/home/ashraya/MD_Project/cluster_3_pca_pdb_classification.csv"
for(k in 1:length(pc.xray$au[,1]))
{
  writestr=paste(k,pc.xray$au[k,1],pc.xray$au[k,2],sep = ",")
  cat(writestr,file = file1name,append = TRUE,sep = "\n")
}

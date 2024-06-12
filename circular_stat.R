files_list<-list.files("/home/ashraya/MD_Project/md_data_from_sneha/1lva_residue_phipsi")
for(j in 1:length(files_list)){
  filename=paste("/home/ashraya/MD_Project/md_data_from_sneha/1lva_residue_phipsi/",files_list[j],sep="")
  res_phipsi<-read.csv(filename,sep = " ",header = FALSE)
  phi_variance<-var.circular(res_phipsi$V3*(pi/180.0))
  psi_variance<-var.circular(res_phipsi$V4*(pi/180.0))
  phi_variance<-(phi_variance*180)/pi
  psi_variance<-(psi_variance*180)/pi
  phi_mean<-mean.circular(res_phipsi$V3*(pi/180.0))
  psi_mean<-mean.circular(res_phipsi$V4*(pi/180.0))
  phi_mean<-(phi_mean*180)/pi
  psi_mean<-(psi_mean*180)/pi
  writestr=paste(files_list[j],phi_variance,psi_variance,phi_mean,psi_mean,sep = ",")
  cat(writestr,file = "/home/ashraya/MD_Project/md_data_from_sneha/1LVA/phipsi_circular_stats_new.csv",append = TRUE,sep = "\n")
}

circ_var<-read.csv("/home/ashraya/MD_Project/md_data_from_kalai/circular_var_nongly.csv")
ggplot(circ_var,aes(x=circ_var$var_psi))+geom_histogram(binwidth = 0.005,colour="black", fill="blue")+labs(title="Variance in Psi angle",x="Circular variance")+scale_x_continuous(breaks = seq(0.007,0.8,0.02))+geom_vline(aes(xintercept=mean(circ_var$var_psi)),colour="red",linetype="dashed")+stat_bin(aes(y=..count.., label=ifelse((..count..)==0,"",..count..)), geom="text", vjust=-1,binwidth=0.005, colour="red")
trans_perresidue<-read.csv("/home/ashraya/MD_Project/md_data_from_sneha/1lva_residue_transitions/53_GLU.csv")
plot(trans_perresidue$frame,trans_perresidue$transition,xaxt='n',yaxt='n',ylim = c(0,5),xlim=c(0,10000),xlab="Frame",ylab="Type of transition",main="GLU429 transitions",pch=16,col = c("dark green", "red","blue","magenta")[as.numeric(trans_perresidue$transition)],xaxs="i",yaxs="i")
xticks<-seq(0,10000,1000)
yticks<-seq(0,5,1)
axis(side=1, at=xticks,labels=T)
axis(side=2, at=yticks,labels=c("","A_TO_A","A_TO_D","D_TO_D","D_TO_A",""))
axis(side=3, at=yticks,labels=F)
axis(side=4, at=yticks,labels=F)

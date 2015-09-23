#documentation start
#=============================================================================
# File data
# creator: Alban Ramette
# acknowledgements: Christiane Hassenrück
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# Used for ARISA quality control
# This function plots the number of bands (OTUs) per PCR grouped by sample to visualize the differences in the number of bands among the replicate PCRs of the same sample.
# input: 
# M - PCR by OTU (peak) table (dataframe), last column containing sample names as factor
# output:
# number, mean and standard deviation of bands per PCR (printed), plot
# dependencies:
# none
#=============================================================================
#documentation end




Plot_Replicate_Variation=function(M) {
  N=as.numeric(droplevels(M[,ncol(M)])) #labels
  U=levels(droplevels(M[,ncol(M)]))
  UL=length(U)
  
  M1=M[,-ncol(M)] #remove the last column of N
  M2=matrix(NA,nrow(M1),ncol(M1))
  
  for(i in 1:nrow(M1)){
    M2[i,]=as.numeric(as.vector(M1[i,]))
  }
  
  M2[M2!=0]<-1
  
  Nbands=apply(M2,1,sum) #how many peaks per sample?
  cat("Nber of bands:",Nbands,"\n")
  cat("Mean= ",mean(Nbands),"; sd=",sd(Nbands),"\n")
  
  
  PeakN=data.frame(cbind(N,Nbands))
  
  # boxplot(Nbands~N,ylim=c(0, max(Nbands)+10),main="Peak number variation per sample")
  
  
  plot(N,Nbands,ylim=c(0, max(Nbands)+10),main="Peak number variation per sample",pch=3,col= 2,xlab="Samples",xaxt="n")
  axis(1,at=seq(1:UL),labels=U,tick=TRUE,hadj=2, cex.axis=0.5, las=2, mgp=c(3, 0.1, 0))
  axis(2,at=seq(from=10,to=max(Nbands),by=10),labels=FALSE,tick=TRUE,padj=0)
  for (i in 1:UL){
    Nb=PeakN[PeakN$N==i,2]
    lines(c(i,i), c(min(Nb),max(Nb)),lty=2)
  }
}
#documentation start
#=============================================================================
# File data
# creator: Alban Ramette
# acknowledgements: modified by Christiane Hassenrück
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# Function to calculate minimum or maximum distance (Bray-Curtis and Jaccard) between replicate PCRs in ARISA
# input: 
# M - PCR by OTU table, last column (G) contains sample names
# fun - function to be used, most commonly max or min
# output:
# list with samples without replication ($no.rep) and table with Bray-Curtis (BC) and Jaccard (JC) distance per sample ($DistRep)
# dependencies:
# requires {vegan}
#=============================================================================
#documentation end

RepDist=function(M,fun){
  U=unique(M[,ncol(M)])
  Nsample=length(U)#how many samples
  result=list(no.rep=c(rep(NA,Nsample)),DistRep=matrix(NA,Nsample,2))
  colnames(result$DistRep)<-c("BC","JC")
  rownames(result$DistRep)=U
  
  for (i in 1:Nsample) {
    Mtemp=M[M[,ncol(M)]==U[i],-ncol(M)]    #Last column = rep #
    if (nrow(Mtemp)==1) { result$no.rep[i]=as.character(U[i])}  #accomodate 1x sample
    if (nrow(Mtemp)!=1){
      DtempBC=vegdist(Mtemp)
      result$DistRep[i,1]=fun(DtempBC)
      DtempJC=vegdist(Mtemp,method="jaccard",binary=TRUE)
      result$DistRep[i,2]=fun(DtempJC)
    }
  }#end of for i
  result$no.rep=result$no.rep[!is.na(result$no.rep)]
  return(result)
}

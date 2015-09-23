#documentation start
#=============================================================================
# File data
# creator: Christiane Hassenrück
# acknowledgements: 
# primary authority: Christiane Hassenrück
# other authorities: 
#=============================================================================
# File contents
# Function to select replicate PCRs that are more dissimilar than a certain cut-off from all other replicate PCRs
# Based on the assumption that if a PCR failed, the distance to all other PCRs of that sample should be high than the cut-off,
#  whereas if a PCR did not fail, the distance to at least one other PCR of that sample should be lower than the cut-off,
#  will only select replicate PCRs that are very dissimilar from the rest - cases where the overall dissimilarity is close to cut-off
#  without a clear pattern of one (or more) outlier PCRs will NOT be selected, manual curation might still be necessary
# input: 
# M - PCR by OTU table, PCR name == rownames(M), last column (G) containing sample names
# cutoff - dissimilarity threshold,  PCRs with a dissimilarity of more than or equal this cut-off are selected
# method - determines whether Bray-Curtis ("BC", default) or Jaccard ("JC") dissimilarity are used
# output:
# vector with names of dissimilar PCRs (rownames(M))
# dependencies:
# requires {vegan}
#=============================================================================
#documentation end


findFailed=function(M, cutoff, method="BC") {
  U=unique(droplevels(M[,ncol(M)]))
  Nsample=length(U)#how many samples
  
  result=vector("list",Nsample)
  
  if(method=="BC"){
    for(i in 1:Nsample) {
      if(length(M$G[droplevels(M$G)==U[i]])!=1){
        Mtemp=M[droplevels(M[,ncol(M)])==U[i],-ncol(M)]
        Dtemp=as.matrix(vegdist(Mtemp), upper=TRUE,diag=TRUE)
        
        for (j in 1:nrow(Mtemp)) {
          if(length(Dtemp[j,][Dtemp[j,]>=cutoff & Dtemp[j,]!=0])==(nrow(Dtemp)-1)){
            result[[i]][j]=rownames(Dtemp)[j]
          }
        }
      }          
    }  
    result.v=unlist(result)
    result.v=result.v[!is.na(result.v)]
    return(result.v)
  }
  
  if(method=="JC"){
    for(i in 1:Nsample) {
      if(length(M$G[droplevels(M$G)==U[i]])!=1){
        Mtemp=M[droplevels(M[,ncol(M)])==U[i],-ncol(M)]
        Dtemp=as.matrix(vegdist(Mtemp, "jaccard", binary=T), upper=TRUE,diag=TRUE)
        
        for (j in 1:nrow(Mtemp)) {
          if(length(Dtemp[j,][Dtemp[j,]>=cutoff & Dtemp[j,]!=0])==(nrow(Dtemp)-1)){
            result[[i]][j]=rownames(Dtemp)[j]
          }
        }
      }    
    }  
    result.v=unlist(result)
    result.v=result.v[!is.na(result.v)]
    return(result.v)
  }
}


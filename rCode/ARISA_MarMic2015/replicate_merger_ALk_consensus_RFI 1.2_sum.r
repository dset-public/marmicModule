#documentation start
#=============================================================================
# File data
# creator: Alban Ramette
# acknowledgements: 
# primary authority:
# other authorities: Christiane Hassenrück
#=============================================================================
# File contents
# function to merge the results of PCR replicates from ARISA into one OTU profile per sample
# input: 
# M - PCR by OTU (peak) table with sample names in the last column (G)
# k - minimum number of times the peak should appear among replicates
# output:
# Sample by OTU table
# dependencies:
# none
#=============================================================================
#documentation end







cat("Replicate merger using the consensus of bands occurring more than 1 per sample\n")
cat("M contains a table of replicates (column rows) x bins (columns) for RFI tables\n")
cat("   => e.g. M=read.table('input.txt',h=TRUE,row.names=1)\n")
cat("The last column contains the sample numbers in increasing numbers\n")
cat("The script was designed to accomodate from 1 to 3 replicates/sample.\n")
cat("k represents the minimum number of times the peak should appear among replicates\n")



Merging_ALk=function(M,k=2){
M=as.data.frame(M)
#ordered list of sample names
U=sort(unique(M[,ncol(M)]))

#how many samples
Nsample=length(U)
Res=matrix(0,Nsample,ncol(M)-1)	#matrix of results
Res=as.data.frame(Res)
rownames(Res)=U
colnames(Res)=colnames(M[,1:(ncol(M)-1)])

for (i in 1:Nsample) {
  Mtemp=M[M[,ncol(M)]==U[i],-ncol(M)]    #Last column = rep #
 # if(is.null(nrow(Mtemp))) {Res[i,]=as.numeric(Mtemp)}	#accomodate 1x sample
  if(nrow(Mtemp)==1) {Res[i,]=as.numeric(Mtemp[1,])}	#accomodate 1x sample
   else(Res[i,]=  apply(Mtemp,2,function(x){ x=as.numeric(as.character(x))
                                        if(length(which(x>0))>=k){sum(x[x!=0])}
                                        else{x=0}
                                      }
             )
   )
}#end of for i
Res=Res[,apply(Res,2,function(x) any(x!=0))]
invisible()
return(Res)
}##################

#e.g.
# Merging_ALk(M,k=2)



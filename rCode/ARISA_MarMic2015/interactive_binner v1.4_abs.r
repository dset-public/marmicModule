# version 1.4
# Copyright (C) 2009 Alban Ramette
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.


interactivebinner=function(D,absolute=F){
cat("interactive binner v.1.4. by A. Ramette\n")

#D<-read.table("input.txt",h=TRUE)

cat("----------------ARISA analyser by A. Ramette----------------\n")
cat('NOTE:   The user needs to set the working directory, as follows:\n')
cat('        for instance:   setwd("c:\\\\R\\\\DIR")\n')
cat('        The result files will be saved in that specific DIR directory\n')
cat('NOTE:   The user needs to import the D table into the R workspace before starting\n')
cat("  	e.g. D=read.table(\"input.txt\",h=TRUE)\n")
cat("  	D is a table with 3 columns:\n")
cat(" 	D[,1] sample name for each band\n")
cat(" 	D[,2] the second consists of band sizes\n")
cat(" 	D[,3] the last consists of area (fluorescence in absolute value\n")
cat("   ***********************************************\n")
cat("	     Variable Window size and Shifting value\n")
cat("------------------------------------------------------------\n")

ANS<-readline("Continue? (y/n)....... ")
if(ANS=="y"){


###################### variable declaration
# 1) range size  to be analyzed
Sm=as.numeric(readline("\nSmallest band size in the range? e.g. 100:\t"))
SM=as.numeric(readline(   "Largest band size in the range? e.g. 1000:\t"))
RFIco=as.numeric(readline("Minimum RFI cutoff value, e.g. 0.09%?\t\t"))
WS=as.numeric(readline( "\nWindow size (bp)? (0.5, 1, 1.5, 2)... \t\t"))
Sh=as.numeric(readline(   "Shift size (bp)? (0.1, 0.5, 1)....... \t\t"))
cat("\n")
PLOT=readline("Plot the results? (y/n)\t\t\t\t")
OUT=readline("Outputting to text files? (y/n)\t\t\t")

#######################
cat("Checking for existence of bands in the defined peak range:",paste(Sm,"-",SM,sep =""),"...\n")

Checkrange=function(D){
  L=length(unique(D[,1]))#how many samples?
  V=as.vector(D[D[,2]<SM&D[,2]>Sm,1]) #samples that are ok
  VD=as.vector(unique(D[,1]))
  Diff=setdiff(VD,V)
    if(length(Diff)!=0) {
	     cat("Size definition problems! The program was stopped.\n",
       "Samples with a problem:\n")
  	   print(Diff)
	     OK<-0
	   }
  if(length(Diff)==0) {
	 cat("OK. No problem of sizes\n")
	 OK<-1
   cat("------------------------------------------------------------\n")
   }
  return(OK)
 }#end of the checkrange ########
OK=Checkrange(D)

if(OK==1){
cat("please wait...\n\n")
#==================
# 1) range size
Ds<-subset(D,D[,2]>=Sm)
Ds<-subset(Ds,Ds[,2]<=SM)
Ds=Ds[order(Ds[,1]),] #get the right order of sample names

#======calculation of the RFI=========================
	# for each sample, subset the area and calculate its sum (Tarea)
Tarea<-tapply(Ds[,3],Ds[,1],sum)
	# get the list of unique names of samples
SampNames<-names(Tarea)


DsRFI=cbind(Ds,RFI=rep(NA,nrow(Ds)))
for (s in 1:length(SampNames)){
  SU=subset(Ds,Ds[,1]==SampNames[s])
  DsRFI[Ds[,1]==SampNames[s],4]=SU[,3]*100/Tarea[s]
}


#keep data based on a RFI cutoff
Dfinal<-subset(DsRFI,DsRFI[,4]>=RFIco)

if(length(unique(Dfinal[,1]))!=length(unique(D[,1]))){
cat(" The script was stopped because the chosen RFI cutoff value has removed too many\n") 
cat("peaks for correct computations. Lower the cutoff value, for instance.\n")}
if(length(unique(Dfinal[,1]))==length(unique(D[,1]))){#RFIco check


#RFI table
TableRFI=xtabs(RFI~Dfinal[,1]+Dfinal[,2],Dfinal)

if(absolute==T){
  TableFI=xtabs(Area~Dfinal[,1]+Dfinal[,2],Dfinal)
}

#01Table
Table01=TableRFI
Table01[Table01!=0]<-1

###################################### binner
#	=> TableRFI is further used for the specific binning

NbF=WS/Sh # Nber of binning Frames


### Create a vector of binning labels

##########################
Shiftpeak=function(Data,Step){####  decrease the column label (peak size) by Step
	Blist=as.numeric(colnames(Data))
	Blist_s=Blist-Step
	Data1=Data
	colnames(Data1)=Blist_s
	return(Data1)                  #return a Table with the same content as the original table but with different labels
}#end of the function shiftpeak
##########################

Binner=function(Table=TableRFI,WS=WS,St=St){
  #for a given TableRFI table, return the binned table
	#conversion of the column names (already shifted) into bins (fixed bin)
	
	 Blist=as.numeric(colnames(Table))
  Convbin=rbind(Peak=Blist,BinNber=rep(NA,length(Blist)),BinStart=rep(NA,length(Blist)))
  #=>3 rows: 1) peak 2) bin Nber  3) Binstart value

  #conversion to bin number for each column
   #E end position vector of the bins
   #S start position vector of the bins
     E=seq(Sm,SM,WS);S=E-WS

     for(j in 1:length(Blist)){
    VAL=Convbin[1,j]
    if(length(which(S<=VAL & VAL<E))==0) {Convbin[2,j]=NA ;Convbin[3,j]=max(E)}
    if(length(which(S<=VAL & VAL<E))!=0) {
        Convbin[2,j]=which(S<=VAL & VAL<E)  #identifies the bin number
        Convbin[3,j]=S[Convbin[2,j]]            #identifies the bin starting number
    }
    }
    #Convbin[2:3,]=round(Convbin[2:3,],0)
    colnames(Table)=Convbin[3,]   #renaming the colnames with BinStart to prepare for binning

   #for each row, add the columns that have the same row1
	U1=as.numeric(unique(colnames(Table))) #vector of unique peak sizes
	U1L=length(U1)
	W1=matrix(NA,nrow(Table),U1L) #matrix of results with binned bands
	rownames(W1)=rownames(Table)
	colnames(W1)=U1

	for (u in 1:U1L){
	Dtemp=Table[,colnames(Table)==U1[u]]
	if(is.null(dim(Dtemp))){W1[,u]=Dtemp}
	else{
	W1[,u]=apply(Dtemp,1,sum) #select the columns with the same peak and sum their RFI per line
	}}

	return(W1) # binned table, wherein RFI are summed by bin
}#end of the function Binner#######################
##############################################



#filling the lists

ShV=seq(0,WS-Sh,Sh) #vector of binning shift values

TableBinVector=vector("list",NbF)	#to produce the list of Tables of shifted peak labels
names(TableBinVector)=ShV  #bin names
for (i in 1:NbF){ #filling the list with the peak relabeled tables
  TableBinVector[[i]]=Shiftpeak(Data=TableRFI,Step=ShV[i])
}


BinnedTables=vector("list",NbF)	#to produce the list of Tables of shifted peak labels
names(BinnedTables)=ShV
for (i in 1:NbF){ #filling the list with the peak relabeled tables
   Tep=   Binner(Table=TableBinVector[[i]],WS=WS,St=Sm-ShV[i]) #St (starts)   100 may be changed 
    colnames(Tep)=as.numeric(colnames(Tep))+Sh*(i-1)    # to correct the names of bins for shifting frame 
  BinnedTables[[i]]= Tep   
  #@need another line here to implement the renaming of bin names as a function of the real bin


}
  

#############################
if(absolute==T){
  
  TableBinVector_abs=vector("list",NbF)  #to produce the list of Tables of shifted peak labels
  names(TableBinVector_abs)=ShV  #bin names
  for (i in 1:NbF){ #filling the list with the peak relabeled tables
    TableBinVector_abs[[i]]=Shiftpeak(Data=TableFI,Step=ShV[i])
  }
  
  
  BinnedTables_abs=vector("list",NbF)  #to produce the list of Tables of shifted peak labels
  names(BinnedTables_abs)=ShV
  for (i in 1:NbF){ #filling the list with the peak relabeled tables
    Tep=   Binner(Table=TableBinVector_abs[[i]],WS=WS,St=Sm-ShV[i]) #St (starts)   100 may be changed 
    colnames(Tep)=as.numeric(colnames(Tep))+Sh*(i-1)    # to correct the names of bins for shifting frame 
    BinnedTables_abs[[i]]= Tep   
    #@need another line here to implement the renaming of bin names as a function of the real bin
    
    
  }
  
}



################################

	#to produce a list of sample correlations
Corlist=lapply(BinnedTables,function(x) as.dist(cor(t(x))))
MeancorVec=sapply(Corlist,mean)




# HELP#####################################################
#     Binning done with a WSize of 1bp with 0.1 bp shifts
#        using the rounding up function
# value: an object "Result" containing the different tables
#        Result[[1]]=    Dfinal      # RFI per peak)
#        Result[[2]]=    TableRFI     # original reformated RFI table)
#        Result[[3]]=    BinnedTables      # List of the binned tables
#        Result[[4]]=    Correlation values
#        Result[[5]]=    Mean correlation per binning frame
#        Result[[6]]=    best frame      # List of the binned tables
#        Result[[7]]=   Summary
# Besides, a plot is produced with a red line indicating the best bin frame (highest mean correlations among samples)
###########################################################

#aggregating the result into one object

Bestbin=names(which(MeancorVec==max(MeancorVec)))
BB=Bestbin

if(length(Bestbin)>1){BB=Bestbin[1]}   #in case several bestbin are found
Result=list(IndivRFI=Dfinal,
            TableRFI=TableRFI,
            BinnedTables=BinnedTables,
            Correlation=Corlist,
            Meancorr=MeancorVec,
            BestFrame=BinnedTables[[BB]],
            BestFrame_abs=BinnedTables_abs[[BB]],
            Summary=rbind(Correlations=round(as.numeric(MeancorVec),2),NberOTUs=round(sapply(BinnedTables,ncol),0))
      )                                                                                                            


if(PLOT=="y"){
  par(mfrow=c(2,1))
  plot(seq(1:NbF),MeancorVec,type="l",ylab="mean inter-sample correlations",xlab="starting point",xaxt="n")
  axis(1,at=seq(1:NbF),ShV)
  #abline(v=as.numeric(Bestbin),col=2)  #no working properly
}

HighOTU=as.numeric(colnames(Result$Summary)[which(as.numeric(Result$Summary[2,])==max(as.numeric(Result$Summary[2,])))])
cat("Best bin frame is: ",Bestbin," (highest mean correlations)\n",sep=" ")

if(PLOT=="y"){
  plot(seq(1:NbF),Result$Summary[2,],type="h",xaxt="n",ylab="Nber OTUs",xlab="Shift",lwd=10,ylim=c(0,max(Result$Summary[2,])))
  axis(1,at=seq(1:NbF),ShV)
 par(mfrow=c(1,1))
}
cat("Max OTU number for frame: ",HighOTU,"\t")
cat("(",max(as.numeric(Result$Summary[2,]))," OTUs)\n",sep="")
cat("(Chosen parameters: WS=",WS,", Shift=",Sh,")\n",sep="")


if(OUT=="y"){
#outputting_______________________________________
  write.table(Dfinal,file="output_D_RFI.txt",append=FALSE,quote=FALSE,row.names =FALSE)
  write.table(t(Table01),file="output_bandlist01.txt",append=FALSE,quote=FALSE)
  write.table(t(TableRFI),file="output_bandlistRFI.txt",append=FALSE,quote=FALSE)
  write.table(t(Result[[6]]),file=paste("output_best.binned.table",Bestbin,".txt",sep=""),append=FALSE,quote=FALSE)
  write.table(t(Result[[7]]),file=paste("output_abs_best.binned.table",Bestbin,".txt",sep=""),append=FALSE,quote=FALSE)
  write.table(Result[[8]],file="output_summary.txt",append=FALSE,quote=FALSE)
}

cat("The results are available in the following object: Result\n\n")
cat("Result summary for each frame:\n")
print(Result[[8]])
return(Result)


cat("(End of calculations)\n")
  }#RFIco OK
 }#ok=1
}# yes

}#end of the function
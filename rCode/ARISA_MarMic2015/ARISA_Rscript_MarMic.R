#input:
# geneMapper output table
# mapping file to match PCR number and sample name
setwd("D:/myTeaching/marmicDsetPublic/trunk/rCode/ARISA_MarMic2015/")

#save.image("MarMic15_ARISA.Rdata")

#reading geneMapper output (export table), in Size column often NAs
Data  <- read.table(
  file = "MarMic15_ARISA_frags.txt",
  h = TRUE, # is there a header?
  fill = TRUE,
  sep = "\t"
  )


#selecting columns of interest
Data1 <- Data[,c("Sample.File.Name","Size","Area")] 

#removing NAs from Size column
Data2 <- Data1[!is.na(Data1$Size),] 

#defining size range (redundant)
Data3 <- Data2[Data2$Size >= 100 & Data2$Size <= 1000,] 

#shortening the names of the samples
ShortNames<-as.character(
  sapply(
    as.character(Data3$Sample.File.Name),
    function(x){
      strsplit(x,"_2015-")[[1]][1]
      }
    )
  )

Data4 <- Data3
Data4$Sample.File.Name <- ShortNames

head(Data4)

write.table(Data4,"for_binning.txt")

#binning
source("interactive_binner v1.4_abs.r")
Binned <- interactivebinner(Data4,absolute=T) #absolute option also outputs the raw peak area values per bin

B <- Binned$BestFrame
B_abs <- Binned$BestFrame_abs
all.equal(rownames(B),rownames(B_abs))

#reading mapping file to match PCR number, i.e. rownames(B), with sample names (SID)
SID0 <- read.table("SID.txt", h=TRUE) 
head(SID0)
PCR <- as.numeric(sapply(as.character(rownames(B)),function(x){strsplit(x,"_")[[1]][2]}))
SID <- SID0[SID0$PCR %in% PCR,]
SID <- SID[match(PCR,SID$PCR),]

#formating input for merging
M <- data.frame(B,G = SID$SID)
M_abs <- data.frame(B_abs,G = SID$SID)

#quality control
require(vegan)

#visual check: negative controls are expected to form halo around samples
NMDS_M <- metaMDS(M[,-ncol(M)])
plot(NMDS_M, display = "sites")

identify(
  NMDS_M$points[,1],
  NMDS_M$points[,2],
  labels=M[,ncol(M)]
  ) #interactive selection of point to print label

#selection of replicate PCRs which are more than a certain cut-off dissimilar from the other replicate PCRs
#Bray_Curtis dissimilarity of <= 0.4 acceptable among replicate PCRs
source("findFailedPCR.R")

failedPCR <- findFailed(M,0.4)

table(
  droplevels(M[failedPCR,"G"])
  )

#visualization of min and max distance between replicate PCRs
source("RepDist.R")
distPCR <- RepDist(M,max)

hist(
  distPCR$DistRep[,1],
  breaks=30,
  main="max dist among rep PCRs"
  ) #BC

hist(
  distPCR$DistRep[,2]
  ) #JC

distPCR <- RepDist(M,min)

hist(distPCR$DistRep[,1],breaks=30) #BC
hist(distPCR$DistRep[,2]) #JC

#variability of number of OTUs in replicate PCRs per sample
source("plotrepvariation_mod.R")
Plot_Replicate_Variation(M)

#removing failed PCRs
M1 <- M[!rownames(M) %in% failedPCR,]
M1_abs <- M_abs[!rownames(M_abs) %in% failedPCR,]

#merging replicate PCRs into 1 sample
#OTU must be present in at least 2 replicate PCRs
source("replicate_merger_ALk_consensus_RFI 1.2.r") #mean of proportions
Merged <- Merging_ALk(M1,k=2)
source("replicate_merger_ALk_consensus_RFI 1.2_sum.r") #sum of raw peak area
Merged_abs <- Merging_ALk(M1_abs,k=2)

write.table(Merged, "Merged_prop.txt", sep="\t", quote = F) #OTUs as proportions/percentages of total peak area
write.table(Merged_abs, "Merged_abs.txt", sep="\t",quote = F) #OTUs as absolute peak area

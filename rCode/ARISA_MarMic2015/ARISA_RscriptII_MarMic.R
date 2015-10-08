#analysis/plotting of ARISA data
#input:
# Merged OTU table from quality control
# table with experimental conditions (time, temperature, spirulina, sample names (SID) identialy to rownames of Merged OTU table)

#reading Merged table
Merged0=read.table("Merged_prop.txt",sep="\t",h=T)
Merged_abs0=read.table("Merged_abs.txt",sep="\t",h=T)

#removing positive control
Merged=Merged0[rownames(Merged0)!="pos",]
Merged_abs=Merged_abs0[rownames(Merged_abs0)!="pos",]

#reading table with experimental conditions
env0=read.table("envData.txt",sep="\t",h=T)
env=env0[match(rownames(Merged),as.character(env0$SID)),]
all.equal(as.character(env$SID),rownames(Merged))
env$int=interaction(env$spirulina,env$temp)
env$col=as.factor(env$temp)
levels(env$col)=c("lightblue","darkgreen","orange","darkred")
env$col=as.character(env$col)
env$col[env$time=="t0"]<-"darkgrey"
env$pch=env$spirulina
levels(env$pch)=c(15,17)
env$pch=as.numeric(as.character(env$pch))
env$label=as.character(env$int)
env$label[env$time=="t0"]<-c("t0")

#seeting plotting parameters
par(xpd=NA,                          # banishes the restriction to only plot within the plotting area. Now also in the margins
    mar=c(5,5,2,10),                  # set margin spaces around the plot c(bottom, left, top, right)
    cex.axis=0.8,
    col.axis="black")            


#alpha diveristy
#number of OTUs
nOTU=rowSums(decostand(Merged,method="pa"))
plot(c(0,max(nOTU))~c(0,max(as.numeric(env$int))),type="n",axes=F,ylab="nOTU",xlab="Experimental condition")
axis(1,at=c(0,unique(as.numeric(env$int))),
     labels=c("t0",levels(env$int)))
axis(2,at=seq(0,max(nOTU),20),labels=T,las=2)
points(nOTU[env$time=="t0"]~rep(0,length(nOTU[env$time=="t0"])),
       pch=env$pch[env$time=="t0"],
       col=env$col[env$time=="t0"],
       cex=1.5)
points(nOTU[env$time!="t0"]~as.numeric(env$int)[env$time!="t0"],
       pch=env$pch[env$time!="t0"],
       col=env$col[env$time!="t0"],
       cex=1.5)

#rarefaction curves
nSEQ=rowSums(Merged_abs)
minSEQ=0
maxSEQ=8000 #adjust value to position when curves level out
plot(maxSEQ,max(nOTU),type="n",xlim=c(minSEQ,10000),ylim=c(0,max(nOTU)))
for(i in 1:nrow(Merged_abs)){
  temp=rarefy(t(Merged_abs)[,i],seq(minSEQ,maxSEQ,100))
  lines(seq(minSEQ,maxSEQ,100),temp,col=env$col[i])
  text(maxSEQ,max(temp),label=rownames(Merged_abs)[i],pos=4,cex=0.7)
}


#beta diversity
#cluster diagram
Merged_clust=hclust(vegdist(Merged))
Merged_clust_pa=hclust(vegdist(Merged,binary=T,method="jaccard"))

plot(Merged_clust,labels=env$label,cex=0.7)
rect.hclust(Merged_clust, h=0.6)
rect.hclust(Merged_clust, h=0.8)
group_0.6=cutree(Merged_clust, h=0.6)
group_0.8=cutree(Merged_clust, h=0.8)

plot(Merged_clust_pa,labels=env$label,cex=0.7)
rect.hclust(Merged_clust_pa, h=0.6)
rect.hclust(Merged_clust_pa, h=0.8)
group_0.6_pa=cutree(Merged_clust_pa, h=0.6)
group_0.8_pa=cutree(Merged_clust_pa, h=0.8)

#NMDS
Merged_nmds=metaMDS(Merged,k=2,trymax=50)
plot(Merged_nmds$points[,1],Merged_nmds$points[,2],type="n",axes=F,ylab="",xlab="")
ordihull(Merged_nmds,groups=group_0.6,lwd=1, lty=1)
ordihull(Merged_nmds,groups=group_0.8,lwd=1, lty=3)
#ordihull(Merged_nmds,groups=group_0.6_pa,lwd=1, lty=1)
#ordihull(Merged_nmds,groups=group_0.8_pa,lwd=1, lty=3)
points(Merged_nmds$points[,1],Merged_nmds$points[,2],pch=env$pch,col=env$col,cex=1.5)

#legend
legend(max(Merged_nmds$points[,1]),max(Merged_nmds$points[,2]),
       legend=c("Time","t0","","Temperature","0°C","5°C","18°C","28°C","","Spirulina","control","fed"),
       col=c("white","darkgrey","white","white","lightblue","darkgreen","orange","darkred","white","white","black","black"),
       pch=c(rep(15,11),17),pt.cex=1.3)





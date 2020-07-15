#### PLS
library(mdatools)
library(pheatmap)
library(ggfortify)

## load chemical data
setwd('/Users/sharp/Desktop/Yasmeen/raw data')
data <- read.csv('chemical-TBMU.csv', row.names = 'Chemical')
group <- read.csv('Chemical grouping.csv')
mdl <- read.csv('mdl.csv',row.names='chemical')
chemData<-data


# filter Data chemicals with many "non-detects"
highDetect<-c()
for(i in 1:nrow(chemData)){
  highDetect<-c(highDetect,sum(chemData[i,]>mdl[i,])/ncol(chemData)>0.75)
}
filterData<-chemData[which(highDetect==1),]

# replace remaining "zero" values with method detection limit
correctedData<-filterData
hasZeros<-apply(filterData, 1, function(x){ any(x==0) })
for(i in names(hasZeros)){
  correctedData[i,correctedData[i,]<mdl[i,]]<-mdl[i,]
}


# Load gene expression data
library(edgeR)
setwd('/Users/sharp/Desktop/Yasmeen/raw data')
filename <- 'PCR-TBMU'
review <- read.csv(paste0(filename,'.csv'), row.names = 'Gene')
review <- as.matrix(review)
for(i in 1:nrow(review)){
  if(length(which(review[i,]==0))){
    review[i,][which(review[i,]==0)] <- mean(review[i,][which(review[i,]!=0)])
  }
}

expData <- 2^(-review)
fcData <-t(apply(normData, 1, function(x) x/mean(x)))
plsGeneData<-fcData


#normalization
sizeFactor <- calcNormFactors(expData, method = "TMM", logratioTrim = 0.5, sumTrim = 0.1) #more strict criteria
normData <- expData
for(i in 1:ncol(normData)) {
  normData[,i] <- expData[,i]/(sum(expData[,i])*sizeFactor[i])
}
plsGeneData<-normData

# Data to use for PLS
finalChem<-log(t(correctedData))
finalGene<-log2(t(plsGeneData[,-c(28)]))
nGenes<-vector()
nComps<-vector()
xVar<-vector()
yVar<-vector()
sMean<-vector()
sMed<-vector()
sMin<-vector()
sMax<-vector()

#pls model-all genes
plsMod<-pls(finalGene, finalChem, center=TRUE, scale=TRUE, cv=1)
nGenes<-c(nGenes,ncol(finalGene))
nComps<-c(nComps,plsMod$ncomp.selected)
#plot top slopes
par(mfrow=c(3,3), mar=c(4,3,1.5,1))
for(i in order(plsMod$cvres$slope[,3],decreasing=TRUE)[1:9]){
  plot(
    plsMod$cvres$y.pred[,3,i]~plsMod$cvres$y.ref[,i],
    xlab=NA,
    ylab=NA,
    main=dimnames(plsMod$calres$y.pred)[[3]][i],
    col='red', pch=16, cex=1.8
  )
  text(plsMod$cvres$y.ref[,i], plsMod$cvres$y.pred[,nComps,i], labels=rownames(finalChem))
  tempLM<-lm(plsMod$cvres$y.pred[,nComps,i]~plsMod$cvres$y.ref[,i])
  abline(tempLM, col="blue", lwd=2)
  abline(0,1, lty=2)
  legend('bottomright', legend = round(plsMod$cvres$slope[i,nComps],2))
}
plsMod.2 <- plsMod

rmseMean<-mean(plsMod.2$cvres$rmse[,nComps[1]])
rmseMed<-median(plsMod.2$cvres$rmse[,nComps[1]])
rmseMin<-min(plsMod.2$cvres$rmse[,nComps[1]])
rmseMax<-max(plsMod.2$cvres$rmse[,nComps[1]])

xVar<-plsMod.2$calres$xdecomp$cumexpvar[nComps[1]]
yVar<-plsMod.2$calres$ydecomp$cumexpvar[nComps[1]]


sMean<-c(sMean, mean(plsMod.2$cvres$slope[,1]))
sMed<-c(sMed, median(plsMod.2$cvres$slope[,1]))
sMin<-c(sMin, min(plsMod.2$cvres$slope[,1]))
sMax<-c(sMax, max(plsMod.2$cvres$slope[,1]))



#pls model-specific VIP gene for each chemical
plsVIPscores<-vipscores(plsMod.2)
hasVIPs<-names(which(apply(plsVIPscores, 2, function(x){any(x>1)})))
#unique Mods based on VIPs per chem
uniqueMods<-list()
uniqueVIPgene<-list()

for(i in hasVIPs){
  try(uniqueVIPgene[[i]]<-finalGene[,plsVIPscores[,i]>1],silent=T)
  try(uniqueMods[[i]]<-pls(uniqueVIPgene[[i]][,], finalChem[,i], center=TRUE, scale=TRUE, cv=1),silent=T)
}
plotOrder<-order(sapply(uniqueMods, function(x){x$cvres$slope[,x$ncomp.selected]}), decreasing = TRUE)[1:12]

nGenes<-c(nGenes,length(hasVIPs))
nComps<-c(nComps,1)

rmseMean<-c(rmseMean, mean(sapply(uniqueMods, function(x) x$cvres$rmse[,nComps[1]])))
rmseMed<-c(rmseMed, median(sapply(uniqueMods, function(x) x$cvres$rmse[,nComps[1]])))
rmseMin<-c(rmseMin, min(sapply(uniqueMods, function(x) x$cvres$rmse[,nComps[1]])))
rmseMax<-c(rmseMax, max(sapply(uniqueMods, function(x) x$cvres$rmse[,nComps[1]])))

tempVar<-median(sapply(uniqueMods, function(x) x$calres$xdecomp$cumexpvar[nComps[1]]))
xVar<-c(xVar, tempVar)

tempVar<-median(sapply(uniqueMods, function(x) x$calres$ydecomp$cumexpvar[nComps[1]]))
yVar<-c(yVar, tempVar)

sMean<-c(sMean, mean(sapply(uniqueMods, function(x) x$cvres$slope[,1])))
sMed<-c(sMed, median(sapply(uniqueMods, function(x) x$cvres$slope[,1])))
sMin<-c(sMin, min(sapply(uniqueMods, function(x) x$cvres$slope[,1])))
sMax<-c(sMax, max(sapply(uniqueMods, function(x) x$cvres$slope[,1])))


par(mfrow=c(4,3), mar=c(4,3,1.5,1))
for(i in plotOrder){
  u_nComps<-uniqueMods[[i]]$ncomp.selected
  dataRange<-c(
    min(c(uniqueMods[[i]]$cvres$y.pred[,u_nComps,1],uniqueMods[[i]]$cvres$y.ref[,1]))-0.2,
    max(c(uniqueMods[[i]]$cvres$y.pred[,u_nComps,1],uniqueMods[[i]]$cvres$y.ref[,1]))+0.2
  )
  plot(
    uniqueMods[[i]]$cvres$y.pred[,u_nComps,1]~uniqueMods[[i]]$cvres$y.ref[,1],
    xlab=NA,
    ylab=NA,
    main=paste0(colnames(finalChem)[i]),
    ylim=dataRange,
    xlim=dataRange,
    col='red', pch=16, cex=1.8
  )
  #slope <- round(uniqueMods[[i]]$cvres$slope[,uniqueMods[[i]]$ncomp.selected],2)
  #mtext(side=1, line=2, cex=0.7, paste0('slope = ', slope, adj=0))
  #r2<-format(round(uniqueMods[[i]]$cvres$r2[,uniqueMods[[i]]$ncomp.selected],2))
  #mtext(side=1, line=2, cex=0.7, bquote(r^2 == .(r2)), adj=1)
  mtext(paste0(
    "slope=", round(uniqueMods[[i]]$cvres$slope[,u_nComps],2),
    " rmse=", round(uniqueMods[[i]]$cvres$rmse[,u_nComps],2)
  ),
  side=1, adj=1, line=2, cex=0.8)
  #text(uniqueMods[[i]]$cvres$y.ref[,1], uniqueMods[[i]]$cvres$y.pred[,u_nComps,1], labels=rownames(finalChem))
  tempLM<-lm(uniqueMods[[i]]$cvres$y.pred[,u_nComps,1]~uniqueMods[[i]]$cvres$y.ref[,1])
  abline(tempLM, col="blue", lwd=2)
  abline(0,1, lty=2)
}


###top genes
VIPchem <- plsVIPscores[,hasVIPs]
topgene <- matrix(nrow=10, ncol=ncol(VIPchem))
colnames(topgene) <- colnames(VIPchem)
write.csv(apply(VIPchem, 2, function(x) names(x[order(x)[length(x):1]][1:10])),
          ' top gene list-TBMU.csv')


### summary of variance
modRes2<-data.frame(
  model=1:2,
  type=c("all genes","VIP gene >1"),
  nGenes=nGenes,
  nComps=nComps,
  gene_Var=xVar,
  chem_Var=yVar,
  rmseMean=rmseMean,
  rmseMed=rmseMed,
  rmseMin=rmseMin,
  rmseMax=rmseMax,
  slopeMean=sMean,
  slopeMed=sMed,
  slopeMin=sMin,
  slopeMax=sMax,
  stringsAsFactors=FALSE
)
modRes2

write.csv(modRes2, 'PLS-TBMU.csv')


## chemical group sum
library(pheatmap)
setwd('/Users/sharp/Desktop/Yasmeen/raw data')
review <- read.csv('chemical-TBMU.csv', row.names = 'Chemical')
group <- read.csv('Chemical grouping.csv')

sum.group <- as.matrix(unique(group[,2]))
sum.mat <- matrix(nrow=11,ncol=ncol(review))
rownames(sum.mat) <- c(sum.group, rownames(review)[1:5])
colnames(sum.mat) <- colnames(review)
for(i in 1:length(sum.group)){
  sum.mat[i, ]<- colSums(review[which(rownames(review) %in% 
                                        group[which(group[,2] %in% sum.group[i]),1]),])
}
sum.mat[7:11,] <- as.matrix(review[1:5,])

### PLS
## Libraries
library(mdatools)
library(edgeR)

## chemical data
singlePAH_Data <-review[-c(1:5),]
sumPAH_Data <- sum.mat[-c(7:11),]
singleMetal_Data <-review[c(1:5),]
sumMetal_Data <- sum.mat[c(7:11),]

## PCR data
library(pheatmap)
library(edgeR)
setwd('/Users/sharp/Desktop/Yasmeen/raw data')
filename <- 'PCR-TBMU-0.1'
review <- read.csv(paste0(filename,'.csv'), row.names = 'Gene')
review <- as.matrix(review)
for(i in 1:nrow(review)){
  if(length(which(review[i,]==0))){
    review[i,][which(review[i,]==0)] <- mean(review[i,][which(review[i,]!=0)])
  }
}
# outlier
a <- which(rownames(review) == 'RPL4')
review[a,30] <- mean(review[a,][which(review[a,]!=0)])
# normalization
expData <- 2^(-review)
sizeFactor <- calcNormFactors(expData, method = "TMM")
normData <- expData
for(i in 1:ncol(normData)) {
  normData[,i] <- expData[,i]/(sum(expData[,i])*sizeFactor[i])
}

qPCRdata <- normData
qPCRdata <- qPCRdata[,which(colnames(qPCRdata) %in% colnames(sumPAH_Data))]
qPCRdata <- t(qPCRdata)

singlePAH_Data <-t(singlePAH_Data)
sumPAH_Data <- t(sumPAH_Data)
singleMetal_Data <-t(singleMetal_Data)
sumMetal_Data <- t(sumMetal_Data)

### Curate Chemical Data 
chemCats<-c("singlePAH", "sumPAH",'singleMetal','sumMetal') # must be spelled the same the .raw data objects

chemData<-list() 
chemPre<-list() 
pre<-vector() # to track filtering results
post<-vector() # to track filtering results
for(i in  chemCats){
  chemData[[i]]<- get(paste(i, "_Data", sep=""))
  rownames(chemData[[i]])<-rownames(chemData[[i]])
  chemPre[[i]]<-chemData[[i]]
  chemData[[i]]<-chemData[[i]][,apply(chemData[[i]], 2, function(x) !all(x==0))] # to return only columns where all vaues are not equal to zero
  chemData[[i]]<-chemData[[i]][,apply(chemData[[i]],2, function(x) sum(x!=0)>3)] # to return only columns chemical is detected in at least 4 of 6 sites
}

### PLS ANALYSIS
# Chemical category to model (choose 1-4 from chemCats) 1="singlePAH", 2="sumPAH", 3="singleMetal", 4="sumMetal"
x <- qPCRdata
y<-chemData[[3]]

plsMods<-list()   # list to hold all models
nGenes<-vector()  # vector to save number of genes used per model
nComps<-vector()  # vector to save number of components used for each model

# MODEL 1: PLS based on all genes
plsMods[[1]]<-pls(x, y, center=TRUE, scale=TRUE, cv=1)
nGenes<-c(nGenes,ncol(qPCRdata))
nComps<-c(nComps,plsMods[[1]]$ncomp.selected)

# MODEL 2: PLS based on all genes, then a VIP test
vipTest<-apply(vipscores(plsMods[[1]]), 1, FUN=function(x) any(x>1))
sum(vipTest)
plsMods[[2]]<-pls(x[,vipTest], y, center=TRUE, scale=TRUE, cv=1)
nGenes<-c(nGenes,sum(vipTest))
nComps<-c(nComps,plsMods[[2]]$ncomp.selected)

### PLS results 
# explained X variance for each model
xVar<-vector()
for(i in 1:length(plsMods)){
  tempVar<-plsMods[[i]]$calres$xdecomp$cumexpvar[nComps[i]]
  xVar<-c(xVar, tempVar)
}

# explained y variance for each model
yVar<-vector()
for(i in 1:length(plsMods)){
  tempVar<-plsMods[[i]]$calres$ydecomp$cumexpvar[nComps[i]]
  yVar<-c(yVar, tempVar)
}

# root mean sqaure error of prediction based on cv results for each model
rmseMean<-vector()
rmseMed<-vector()
rmseMin<-vector()
rmseMax<-vector()
for(i in 1:length(plsMods)){
  rmseMean<-c(rmseMean, mean(plsMods[[i]]$cvres$rmse[,nComps[i]]))
  rmseMed<-c(rmseMed, median(plsMods[[i]]$cvres$rmse[,nComps[i]]))
  rmseMin<-c(rmseMin, min(plsMods[[i]]$cvres$rmse[,nComps[i]]))
  rmseMax<-c(rmseMax, max(plsMods[[i]]$cvres$rmse[,nComps[i]]))
}

# slope of prediction based on cv results for each model
sMean<-vector()
sMed<-vector()
sMin<-vector()
sMax<-vector()
for(i in 1:length(plsMods)){
  sMean<-c(sMean, mean(plsMods[[i]]$cvres$slope[,3]))
  sMed<-c(sMed, median(plsMods[[i]]$cvres$slope[,3]))
  sMin<-c(sMin, min(plsMods[[i]]$cvres$slope[,3]))
  sMax<-c(sMax, max(plsMods[[i]]$cvres$slope[,3]))
}


modRes2<-data.frame(
  model=1:2,
  type=c("all genes","all VIP>1"),
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



### Explore and plot results
bestMod<-1

# chemical that has lowest RMSE
which(plsMods[[bestMod]]$cvres$rmse[,3]==min(plsMods[[bestMod]]$cvres$rmse[,3]))

# chemcail that has best prediction slope
which(plsMods[[bestMod]]$cvres$slope[,3]==max(plsMods[[bestMod]]$cvres$slope[,3]))

# plot the top 9 RMSEs
par(mfrow=c(3,3), mar=c(4,3,1.5,1))
for(i in order(plsMods[[bestMod]]$cvres$rmse[,3])[1:9]){
  plot(
    plsMods[[bestMod]]$cvres$y.pred[,3,i]~plsMods[[bestMod]]$cvres$y.ref[,i],
    xlab=NA,
    ylab=NA,
    main=dimnames(plsMods[[bestMod]]$calres$y.pred)[[3]][i],
    col='red', pch=16, cex=1.8)
  tempLM<-lm(plsMods[[bestMod]]$cvres$y.pred[,3,i]~plsMods[[bestMod]]$cvres$y.ref[,i])
  abline(tempLM, col="blue", lwd=2)
  abline(0,1, lty=2)
}


# plot the top 9 slopes
par(mfrow=c(2,3), mar=c(4,3,1.5,1))
for(i in order(plsMods[[bestMod]]$cvres$slope[,3], decreasing=TRUE)[1:6]){
  plot(
    plsMods[[bestMod]]$cvres$y.pred[,3,i]~plsMods[[bestMod]]$cvres$y.ref[,i],
    xlab=NA,
    ylab=NA,
    main=dimnames(plsMods[[bestMod]]$calres$y.pred)[[3]][i],
    col='red', pch=16, cex=1.8, cex.main = 1)
  tempLM<-lm(plsMods[[bestMod]]$cvres$y.pred[,3,i]~plsMods[[bestMod]]$cvres$y.ref[,i])
  abline(tempLM, col="blue", lwd=2)
  abline(0,1, lty=2)
  legend('bottomright', legend = round(plsMods[[bestMod]]$cvres$slope[i,3],2))
}

top.gene <- matrix(nrow=10,ncol=9)
m1 <- order(plsMods[[bestMod]]$cvres$slope[,3], decreasing=TRUE)[1:9]
colnames(top.gene) <- dimnames(plsMods[[bestMod]]$coeffs$values)[[3]][m1]
for(i in 1:9){
  m2 <- sort(abs(plsMods[[bestMod]]$coeffs$values[,,m1[i]])[,1])
  top.gene[,i] <- as.matrix(names(m2[length(m2):(length(m2)-9)]))
}
top.gene



## chemical group sum
library(pheatmap)
library(ggfortify)

setwd('/Users/sharp/Desktop/Yasmeen/raw data')
data <- read.csv('chemical-TBMU.csv', row.names = 'Chemical')
group <- read.csv('Chemical grouping.csv')
mdl <- read.csv('mdl.csv',row.names='chemical')
chemData<-data


# filter Data chemicals with many "non-detects"
highDetect<-c()
for(i in rownames(chemData)){
  highDetect<-c(highDetect,sum(chemData[i,]>mdl[i,])/ncol(chemData)>0.25)
}
filterData<-chemData[which(highDetect==1),]

# replace remaining "zero" values with method detection limit
correctedData<-filterData
hasZeros<-apply(filterData, 1, function(x){ any(x==0) })
for(i in names(hasZeros)){
  correctedData[i,correctedData[i,]<mdl[i,]]<-mdl[i,]
}

# PCA
pcaData <- t(log10(correctedData))
pcaData <- t(correctedData)
rownames(pcaData)<-sapply(strsplit(rownames(pcaData),'U.'),function(x)x[2])

pcaMod<-prcomp(pcaData)

PCA<-pcaMod
top5_pc1<-names(PCA$rotation[order(abs(PCA$rotation[,1]), decreasing=TRUE),1])[1:5]
top5_pc2<-names(PCA$rotation[order(abs(PCA$rotation[,2]), decreasing=TRUE),2])[1:5]
top5<-unique(c(top5_pc1, top5_pc2))

PCA2<-PCA
PCA2$rotation<-PCA2$rotation[top5,]
PCA2$center<-PCA2$center[top5]
PCA2$scale<-PCA2$scale[top5]

p<-autoplot(PCA2, data=pcaData,label=F, 
            shape = 16, size = 4,  colour = 'red',
            loadings = TRUE, loadings.label = TRUE,
            loadings.colour = 'darkgrey',
            loadings.label.colour = 'darkgrey',
            loadings.label.size = 4,xlim=c(-0.5,0.5),
            main='PCA of TBMU chemical burdens') + 
  geom_text(aes(label = rownames(pcaData)), size = 4,
            nudge_y = rep(0.025,dim(pcaData)[1])) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA)) + 
  geom_hline(aes(yintercept=0), colour="black", linetype="dashed") +  
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed")
p

write.csv(cbind(top5_pc1,top5_pc2),'top_PCA of TBMU chemical burdens.csv')


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

#normalization
sizeFactor <- calcNormFactors(expData, method = "TMM", logratioTrim = 0.5, sumTrim = 0.1) #more strict criteria
normData <- expData
for(i in 1:ncol(normData)) {
  normData[,i] <- expData[,i]/(sum(expData[,i])*sizeFactor[i])
}
fcData <-t(apply(normData, 1, function(x) x/mean(x)))
plsGeneData<-log(fcData,2)


# PCA
pcaData <- t(plsGeneData)
rownames(pcaData)<-sapply(strsplit(rownames(pcaData),'U.'),function(x)x[2])
pcaMod<-prcomp(pcaData[,-ncol(pcaData)],
               center=TRUE,
               scale=TRUE)

PCA<-pcaMod
top5_pc1<-names(PCA$rotation[order(abs(PCA$rotation[,1]), decreasing=TRUE),1])[1:5]
top5_pc2<-names(PCA$rotation[order(abs(PCA$rotation[,2]), decreasing=TRUE),2])[1:5]
top5<-unique(c(top5_pc1, top5_pc2))

PCA2<-PCA
PCA2$rotation<-PCA2$rotation[top5,]
PCA2$center<-PCA2$center[top5]
PCA2$scale<-PCA2$scale[top5]

p<-autoplot(PCA2, data=pcaData,label=F, 
            shape = 16, size = 4,  colour = 'red',
            loadings = TRUE, loadings.label = TRUE,
            loadings.colour = 'darkgrey',
            loadings.label.colour = 'darkgrey',
            loadings.label.size = 4,xlim=c(-0.5,0.5),
            main='PCA of TBMU gene expressions') + 
  geom_text(aes(label = rownames(pcaData)), size = 4,
            nudge_y = rep(0.025,dim(pcaData)[1])) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA)) + 
  geom_hline(aes(yintercept=0), colour="black", linetype="dashed") +  
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed")
p


write.csv(cbind(top5_pc1,top5_pc2),'top_PCA of TBMU gene expressions.csv')





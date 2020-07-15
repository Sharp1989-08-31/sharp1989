library(pheatmap)
library(edgeR)
library(ggplot2)
setwd('/Users/sharp/Desktop/Yasmeen/raw data')
filename <- 'PCR-TBMU'
review <- read.csv(paste0(filename,'.csv'), row.names = 'Gene')
review <- as.matrix(review)
for(i in 1:nrow(review)){
  if(length(which(review[i,]==0))){
    review[i,][which(review[i,]==0)] <- mean(review[i,][which(review[i,]!=0)])
  }
}

#outlier
a <- which(rownames(review) == 'RPL4');review[a,30] <- mean(review[a,][which(review[a,]!=0)])


expData <- 2^(-review)

sizeFactor <- calcNormFactors(expData, method = "TMM")
normData <- expData
for(i in 1:ncol(normData)) {
  normData[,i] <- expData[,i]/(sum(expData[,i])*sizeFactor[i])
}

fcData <-log(t(apply(normData, 1, function(x) x/mean(x))),2)


###Heatmap
library(circlize)
library(ComplexHeatmap)

colnames(fcData)<-sapply(strsplit(colnames(fcData),'U.'),function(x)x[2])

color_fun<- colorRamp2(c(min(fcData),mean(fcData),max(fcData)), c("green", "black", "red"))

dcols = dist(t(fcData), method = "manhattan")
drows = dist(fcData, method = "manhattan")

ht1<-Heatmap(as.matrix(fcData),name = "log2 gene expression",
             show_row_dend=F,
             clustering_distance_columns = dcols,
             clustering_distance_rows = drows,
             column_title = 'TBMU Chemical Residues',
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 8),
             #rect_gp = gpar(col = "black", lwd = 1),
             col = color_fun)
ht1

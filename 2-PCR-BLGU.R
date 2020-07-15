library(pheatmap)
library(edgeR)
setwd('/Users/sharp/Desktop/Yasmeen/raw data')
filename <- 'PCR-BLGU'
review <- read.csv(paste0(filename,'.csv'), row.names = 'Gene')
review <- as.matrix(review)
for(i in 1:nrow(review)){
  if(length(which(review[i,]==0))){
    review[i,][which(review[i,]==0)] <- mean(review[i,][which(review[i,]!=0)])
  }
}

#outlier
#a <- which(rownames(review) == 'ALDH1A1')
#review[a,27] <- mean(review[a,][which(review[a,]!=0)])

expData <- 2^(-review)
sizeFactor <- calcNormFactors(expData, method = "TMM")
normData <- expData
for(i in 1:ncol(normData)) {
  normData[,i] <- expData[,i]/(sum(expData[,i])*sizeFactor[i])
}
#boxplot(-log2(normData))
fcData <-log(t(apply(normData, 1, function(x) x/mean(x))),2)



###Heatmap
library(circlize)
library(ComplexHeatmap)
ann_colors = list(Site=c(Site_1='cornflowerblue',Site_2='yellow'),
                  Categroy=c('Trace Elements'='grey',
                             'Parent PACs'='black',
                             'Alkylated PACs'='white'))

colnames(fcData)<-sapply(strsplit(colnames(fcData),'U.'),function(x)x[2])
Site<-as.matrix(c(rep('Site_1',16),rep('Site_2', 14)))
rownames(Site)<-colnames(fcData)
annotation_col <- data.frame(
  Site = Site
)

ha_col<-HeatmapAnnotation(Site=Site,
                          col=ann_colors[1],border = TRUE,
                          gp = gpar(col = "black"))


color_fun<- colorRamp2(c(min(fcData),mean(fcData),max(fcData)), c("green", "black", "red"))

dcols = dist(t(fcData), method = "manhattan")
drows = dist(fcData, method = "manhattan")

ht1<-Heatmap(as.matrix(fcData),name = "log2 gene expression",
             top_annotation = ha_col,
             show_row_dend=F,
             clustering_distance_columns = dcols,
             clustering_distance_rows = drows,
             column_title = 'BLGU Chemical Residues',
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 8),
             #rect_gp = gpar(col = "black", lwd = 1),
             col = color_fun)
ht1

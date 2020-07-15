## chemical group sum
library(pheatmap)
library(ggfortify)
library(ComplexHeatmap)
library(circlize)

setwd('/Users/sharp/Desktop/Yasmeen/raw data')
data <- read.csv('chemical-BLGU.csv', row.names = 'Chemical')
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
colSums(correctedData)
correctedData <- log10(correctedData)


###Heatmap
ann_colors = list(Site=c(Site_1='darkgreen',Site_2='yellow'),
                  Group=c(PLMW='pink',
                          PHMW='red',
                          PHET='purple',
                          ALMW='lightblue',
                          AHMW='blue',
                          AHET='darkblue',
                          'Trace Elements'='lightgreen'),
                  Categroy=c('Trace Elements'='grey',
                             'Parent PACs'='black',
                             'Alkylated PACs'='white'))

colnames(correctedData)<-sapply(strsplit(colnames(correctedData),'U.'),function(x)x[2])
Site<-as.matrix(c(rep('Site_1',17),rep('Site_2', 11)))
rownames(Site)<-colnames(correctedData)
annotation_col <- data.frame(
  Site = Site
)

ha_col<-HeatmapAnnotation(Site=Site,
                          col=ann_colors[1],border = TRUE)

setwd('/Users/sharp/Desktop/Yasmeen/raw data')
group<-as.matrix(read.csv('Chemical grouping.csv'))
group<-group[which(group[,1]%in%rownames(as.matrix(correctedData))),]
rownames(group)<-group[,1]
group<-group[order(group[,2])[nrow(group):1],]
Group<-as.matrix(group[,2])

Categroy<-as.matrix(group[,3])

ha_row<-rowAnnotation(Group=Group,
                      Categroy=Categroy,
                      col=ann_colors[2:3],
                      border = TRUE)

fcData<-as.matrix(correctedData)
color_fun<- colorRamp2(c(min(fcData),mean(fcData),max(fcData)), c("cornflowerblue", "#EEEEEE", "red"))

dcols = dist(t(correctedData), method = "manhattan")

row_order<-c()
for(i in 1:nrow(correctedData)){
  row_order<-c(row_order,which(rownames(correctedData)%in%rownames(group)[i]))
}
correctedData<-correctedData[row_order,]
ht1<-Heatmap(as.matrix(correctedData),name = "log10(ng/g)",
             top_annotation = ha_col,
             right_annotation = ha_row,
             show_row_dend=F,
             clustering_distance_columns = dcols,
             column_title = 'BLGU Chemical Residues',
             cluster_rows = F,
             column_names_gp = gpar(fontsize = 10),
             row_names_gp = gpar(fontsize = 8),
             col = color_fun)
ht1



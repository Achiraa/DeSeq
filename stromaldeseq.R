setwd("C:/geo_data/stromalDeseq")
R.version

#installing required libraries
library(DESeq2)
library(tidyverse)
library(BiocGenerics)
library(S4Vectors)
library(stats4)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(gplots)
library(ComplexHeatmap)
library(circlize)

#reading the sample info of the experiment containing samplenames, cell info
data1=read.csv("sampleinfo1.csv")
dim(data1)
data1
head(data1)
View(data1)

#reading the csv file with gene counts and ensembleID
countdat=read.csv("finalcsvfile1.csv")
countdat1=countdat[,-1]
countdat1
row.names(countdat1)=countdat[,1]
head(countdat1)

dim(countdat1)
View(countdat1)

#checking whether the rownames and column names in both the files are similar
all(colnames(countdat1) %in% rownames(data1))
all(colnames(countdat1) == rownames(data1))

#converting the count data into different forms like class, assays etc
#used to store the input values, intermediate calculations and results of an analysis of differential expression
dds=DESeqDataSetFromMatrix(countData = countdat1,colData = data1, design = ~ Cell.type)
dds
View(dds)
View(counts(dds))

#only keeping the genes which have at least 10 counts
keep= rowSums(counts(dds)) >= 10
keep
dds1=dds[keep,]
dds1

#Generate the normalized counts, it uses he median of ratios method of normalization
ddsnew=estimateSizeFactors(dds1)
ddsnew

#checking the normalized counts
sizeFactors(ddsnew)

#saving the normalized counts
#retrieve  the  normalized  counts  matrix  from  dds,  we  use  the  counts()  function  and  
#add  the  argument  normalized=TRUE.
nm1=counts(ddsnew,normalized=TRUE)
nm1
wt1=write.table(nm1, file=" normalized_countsstromal.txt", sep="\t", quote=F,col.names=NA)
wt1

#To  improve  the  distances/clustering  for  the  PCA  and  heirarchical  clustering  visualization  methods,  
#we  need  to  moderate the variance across the mean by applying the rlog transformation to the normalized counts.
rld=rlog(ddsnew, blind=TRUE)
rld

#Differential expression analysis with DESeq2
newdds=DESeq(ddsnew)
newdds

ddsnew$Treatment=relevel(ddsnew$Treatment,ref = "Untreated")
ddsnew$Treatment

#gives us the different column i.e. log2foldchange, baseMean, Pval, padj
resultss=results(newdds)
resultss

#gives us the info regarding how many genes are up and down regulated
summary(resultss)
saveRDS(resultss,file = "resultss.rds")
xy=write.csv(resultss,file = "resuts1.csv",row.names = TRUE,append = TRUE)
xy

#taking FDR>0.01 or psdj>0.01
result.01=results(newdds, alpha=0.01)
result.01
summary(result.01)

#getting the most differentiated expressed genes and saving in CSV
resultsort=result.01[order(result.01$padj),]
topdeseq=resultsort[1:8150,]
write.csv(topdeseq,file = "topdeseqgenes8150.csv")


#PCA plot
plotPCA(rld, intgroup="Cell.type")
plotPCA(rld, intgroup="Treatment")

#MA Plot
plotMA(resultss)
plotMA(result.01)
#idx=identify(resultss$baseMean, resultss$log2FoldChange)
#rownames(resultss)[idx]

#volcano
v=EnhancedVolcano(result.01,lab=rownames(newdds),x="log2FoldChange",y="padj")
v

#heatmap
topgenes=head(rownames(resultsort),8150)
mat1=assay(newdds)[topgenes,]
mat1=mat1 -rowMeans(mat1)
heatmap(mat1)









































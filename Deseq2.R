#setting up the working directories
setwd("C:/geo_data/Deseq")

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
data=read.csv("sampleinfo.csv")
dim(data)
data
data1=as.matrix(data)
head(data1)
View(data1)

#reading the csv file with gene counts and ensembleID
countdata=read.csv("finalcsvfile.csv")
countdata1=countdata[,-1]
countdata1
row.names(countdata1)=countdata[,1]
countnew=as.matrix(countdata1)
head(countnew)
dim(countnew)
View(countnew)

#checking whether the rownames and column names in both the files are similar
all(colnames(countnew) %in% rownames(data1))
all(colnames(countnew) == rownames(data1))

#replacing NA with 0
countnew[is.na(countnew)]=0

#converting the count data into different forms like class, assays etc
#used to store the input values, intermediate calculations and results of an analysis of differential expressionds
ds=DESeqDataSetFromMatrix(countData = countnew,colData = data1,design = ~Treatment)
ds
View(ds)
View(counts(ds))

#only keeping the genes which have at least 10 counts
keep= rowSums(counts(ds)) >= 10
keep
ds=ds[keep,]
ds

#Generate the normalized counts, it uses he median of ratios method of normalization
ds=estimateSizeFactors(ds)
ds

#checking the normalized counts
sizeFactors(ds)

#saving the normalized counts
#retrieve  the  normalized  counts  matrix  from  dds,  we  use  the  counts()  function  and  
#add  the  argument  normalized=TRUE.
nm=counts(ds,normalized=TRUE)
nm
wt=write.table(nm, file=" normalized_counts.txt", sep="\t", quote=F,col.names=NA)
wt

#To  improve  the  distances/clustering  for  the  PCA  and  heirarchical  clustering  visualization  methods,  
#we  need  to  moderate the variance across the mean by applying the rlog transformation to the normalized counts.
rld=rlog(ds, blind=TRUE)
rld

#Differential expression analysis with DESeq2
newds=DESeq(ds)
newds

ds$Treatment=relevel(ds$Treatment,ref = "Untreated")
ds$Treatment

#gives us the different column i.e. log2foldchange, baseMean, Pval, padj
resu=results(newds)
resu

#gives us the info regarding how many genes are up and down regulated
summary(resu)
saveRDS(resu,file = "resu.rds")
xy=write.csv(resu,file = "resut1.csv",row.names = TRUE,append = TRUE)
xy

#taking FDR>0.01 or psdj>0.01
res1=results(newds,contrast = list("Treatment_Treated_vs_Untreated"))
res1
resu.01=results(newds, alpha=0.01)
resu.01
summary(resu.01)

#getting the most differentiated expressed genes and saving in CSV
resusort=resu[order(resu$padj),]
topdeseq=resusort[1:843,]
write.csv(topdeseq,file = "topdeseqgenes.csv")


#PCA plot
plotPCA(rld, intgroup="Cell.type")
plotPCA(rld, intgroup="Treatment")

#MA Plot
plotMA(resu)
plotMA(resu.01)
idx=identify(resu$baseMean, resu$log2FoldChange)
rownames(resu)[idx]

#volcano
v=EnhancedVolcano(resu,lab=rownames(newds),x="log2FoldChange",y="padj")
v

#heatmap
topgenes=head(rownames(resusort),843)
mat1=assay(ds)[topgenes,]
mat1=mat1 -rowMeans(mat1)
heatmap(mat1)




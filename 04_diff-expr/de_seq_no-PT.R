setwd("~/grid/shelf/genomics/projects/ntabacum/analysis/06_prepare_data_for_TBro/02_DESeq")

source("https://bioconductor.org/biocLite.R")
biocLite("tximport")
install.packages("readr")

library(tximport)
library(readr)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(ggplot2)

#reading salmon output
sampleNames=c("OV-A","OV-B","OV-C","OV-D","OV-E",
              "OVF-A","OVF-B","OVF-C",
              "OVP-A","OVP-B","OVP-C")
files=paste("./quants/",sampleNames,".quant.sf",sep="")
names(files)=sampleNames
files
file.exists(files)
tx2gene=read.table("transcript2unigene.IDs", head=F)
colnames(tx2gene)=c("TXNAME","GENEID")
head(tx2gene)

txi.salmon=tximport(files, type="salmon", tx2gene=tx2gene, reader=read_tsv)
head(txi.salmon$counts)

#moving to DESeq
coldata=data.frame(colnames(txi.salmon$counts))
colnames(coldata)=c("sample")
coldata$tissue=c("OV","OV","OV","OV","OV",
                  "OVF","OVF","OVF",
                  "OVP","OVP","OVP")

head(coldata)


dds=DESeqDataSetFromTximport(txi.salmon, coldata, ~tissue)

#exploratory analysis and visualization
# deleting rows with only 0 counts
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

#rlog transformation
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)

## par( mfrow = c( 1, 2 ) )
## dds <- estimateSizeFactors(dds)
## plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
##           pch=16, cex=0.3)
## plot(assay(rld)[,1:2],
##           pch=16, cex=0.3)


## #print normalized counts
## write.table(file="normalized_counts.tab", counts(dds, normalized=T),quote=F, sep="\t")

#sample distances
sampleDists <- dist( t( assay(rld) ) )
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- sampleNames
colnames(sampleDistMatrix) <- NULL
sampleDistMatrix
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
                  clustering_distance_rows=sampleDists,
                  clustering_distance_cols=sampleDists,
                  col=colors)

pdf("plot-PCA_no-PT.pdf")
plotPCA(rld, intgroup = c("tissue"))
dev.off()

## pcaData <- plotPCA(rld, intgroup = c( "sample","tissue"), returnData=TRUE)
## percentVar <- round(100 * attr(pcaData, "percentVar"))
## #pdf("PCA.pdf")
## ggplot(pcaData, aes(PC1, PC2, label=sample, colour=tissue)) +
##     xlab(paste0("PC1: ",percentVar[1],"% variance")) +
##     ylab(paste0("PC2: ",percentVar[2],"% variance")) +
##     ylim(-30,30)+geom_text(size=3)
## #dev.off()



#Differential expression analysis
dds <- DESeq(dds)

## (res <- results(dds))
## summary(res)
## res

#geat base means for each tissue
baseMeanPerLvl <- sapply( levels(dds$tissue), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$tissue == lvl] ) )
head(baseMeanPerLvl)


resPT_OV <- results(dds, contrast=c("tissue","PT","OV"))
summary(resPT_OV)
resPT_OV$baseMeanA=baseMeanPerLvl[,4]
resPT_OV$baseMeanB=baseMeanPerLvl[,1]
head(resPT_OV)
(resPT_OVP <- results(dds, contrast=c("tissue","PT","OVP")))
summary(resPT_OVP)
resPT_OVP$baseMeanA=baseMeanPerLvl[,4]
resPT_OVP$baseMeanB=baseMeanPerLvl[,3]
head(resPT_OVP)
(resPT_OVF <- results(dds, contrast=c("tissue","PT","OVF")))
summary(resPT_OVF)
resPT_OVF$baseMeanA=baseMeanPerLvl[,4]
resPT_OVF$baseMeanB=baseMeanPerLvl[,2]
(resOV_OVF <- results(dds, contrast=c("tissue","OV","OVF")))
summary(resOV_OVF)
resOV_OVF$baseMeanA=baseMeanPerLvl[,1]
resOV_OVF$baseMeanB=baseMeanPerLvl[,2]
(resOV_OVP <- results(dds, contrast=c("tissue","OV","OVP")))
summary(resOV_OVP)
resOV_OVP$baseMeanA=baseMeanPerLvl[,1]
resOV_OVP$baseMeanB=baseMeanPerLvl[,3]
(resOVP_OVF <- results(dds, contrast=c("tissue","OVP","OVF")))
summary(resOVP_OVF)
resOVP_OVF$baseMeanA=baseMeanPerLvl[,3]
resOVP_OVF$baseMeanB=baseMeanPerLvl[,2]






## res.05 <- results(dds, alpha=.05)
## table(res.05$padj < .05)
## resLFC1 <- results(dds, lfcThreshold=1)
## table(resLFC1$padj < 0.1)

## plotMA(res, ylim=c(-7,15))

## topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)
## mat <- assay(rld)[ topVarGenes, ]
## mat <- mat - rowMeans(mat)
## df <- as.data.frame(colData(rld)[,c("tissue","condition")])
## pheatmap(mat, annotation_col=df)


#exporting results
## resOrdered <- res[order(res$padj),]
## head(resOrdered)
## resOrderedDF <- as.data.frame(resOrdered)
## write.csv(resOrderedDF, file="results.csv")


resOrderedPT_OV <- resPT_OV[order(resPT_OV$padj),]
head(resOrderedPT_OV)
resOrderedDFPT_OV <- as.data.frame(resOrderedPT_OV)
write.csv(resOrderedDFPT_OV, file="results_PT-OV.csv")

resOrderedPT_OVF <- resPT_OVF[order(resPT_OVF$padj),]
head(resOrderedPT_OVF)
resOrderedDFPT_OVF <- as.data.frame(resOrderedPT_OVF)
write.csv(resOrderedDFPT_OVF, file="results_PT-OVF.csv")

resOrderedPT_OVP <- resPT_OVP[order(resPT_OVP$padj),]
head(resOrderedPT_OVP)
resOrderedDFPT_OVP <- as.data.frame(resOrderedPT_OVP)
write.csv(resOrderedDFPT_OVP, file="results_PT-OVP.csv")

resOrderedOV_OVP <- resOV_OVP[order(resOV_OVP$padj),]
head(resOrderedOV_OVP)
resOrderedDFOV_OVP <- as.data.frame(resOrderedOV_OVP)
write.csv(resOrderedDFOV_OVP, file="results_OV-OVP.csv")

resOrderedOV_OVF <- resOV_OVF[order(resOV_OVF$padj),]
head(resOrderedOV_OVF)
resOrderedDFOV_OVF <- as.data.frame(resOrderedOV_OVF)
write.csv(resOrderedDFOV_OVF, file="results_OV-OVF.csv")

resOrderedOVP_OVF <- resOVP_OVF[order(resOVP_OVF$padj),]
head(resOrderedOVP_OVF)
resOrderedDFOVP_OVF <- as.data.frame(resOrderedOVP_OVF)
write.csv(resOrderedDFOVP_OVF, file="results_OVP-OVF.csv")



# actual code
library(RColorBrewer)
library(gplots)
library(DESeq2, quietly = TRUE)
setwd("/Users/barrantes/anshu/temp")
exprTable01 <- read.csv("set01_expression_counts.tsv", header=TRUE, sep="\t")
exprTable02 <- read.csv("set02_expression_counts.tsv", header=TRUE, sep="\t")
countsSet01A <- exprTable01[,c(3:4,7:8,9:10)]
countsSet02B <- exprTable02[,c(4:5,7:8,9:10)]
countsSet02C <- exprTable02[,c(4:5,7:8,10:11)]
rownames(countsSet01A) <- exprTable01[,c("target_id")]
rownames(countsSet02B) <- exprTable02[,c("target_id")]
rownames(countsSet02C) <- exprTable02[,c("target_id")]
minCount <- 2700

# set01, group A
countsAboveThr <- subset(countsSet01A, rowSums(countsSet01A[,c(1:6)]) > minCount)
samples <- data.frame(row.names=c("A1", "A2", "C2", "C3", "P1", "P2"), condition=as.factor(c(rep("A",2),rep("C",2),rep("P",2))))
bckCDS2 <- DESeqDataSetFromMatrix(countData = countsAboveThr, colData=samples, design=~condition)
bckCDS2_1 <- DESeq(bckCDS2)
rld <- rlog(bckCDS2_1, blind=FALSE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("set01.distanceheatmap21.20170308.pdf")
# density.info="none",  # turns off density plot inside color legend
# trace="none",         # turns off trace lines inside the heat map
# dendrogram="row",     # only draw a row dendrogram
heatmap.2(sampleDistMatrix, col = rev(colors), main = "Set 01 sample distances", density.info="none", trace="none", dendrogram="row")
dev.off()

# set02, group B
countsAboveThr <- subset(countsSet02B, rowSums(countsSet02B[,c(1:6)]) > minCount)
samples <- data.frame(row.names=c("A2", "A3", "C2", "C3", "P1", "P2"), condition=as.factor(c(rep("A",2),rep("C",2),rep("P",2))))
bckCDS2 <- DESeqDataSetFromMatrix(countData = countsAboveThr, colData=samples, design=~condition)
bckCDS2_1 <- DESeq(bckCDS2)
rld <- rlog(bckCDS2_1, blind=FALSE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("set02.distanceheatmap22.20170308.pdf")
heatmap.2(sampleDistMatrix, col = rev(colors), main = "Set 02 sample distances #1", density.info="none", trace="none", dendrogram="row")
dev.off()

# set02, group C
countsAboveThr <- subset(countsSet02C, rowSums(countsSet02C[,c(1:6)]) > minCount)
samples <- data.frame(row.names=c("A2", "A3", "C2", "C3", "P2", "P3"), condition=as.factor(c(rep("A",2),rep("C",2),rep("P",2))))
bckCDS2 <- DESeqDataSetFromMatrix(countData = countsAboveThr, colData=samples, design=~condition)
bckCDS2_1 <- DESeq(bckCDS2)
rld <- rlog(bckCDS2_1, blind=FALSE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf("set02.distanceheatmap23.20170308.pdf")
heatmap.2(sampleDistMatrix, col = rev(colors), main = "Set 02 sample distances #2", density.info="none", trace="none", dendrogram="row")
dev.off()




library(DESeq2, quietly = TRUE)
library(xlsx) # saving multiple tables into the same Excel Workbook
library(RColorBrewer) # for brewer.pal
library(gplots) # for heatmap.2

setwd("/Users/barrantes/anshu/temp")
exprTable01 <- read.csv("set01_expression_counts.tsv", header=TRUE, sep="\t")
exprTable02 <- read.csv("set02_expression_counts.tsv", header=TRUE, sep="\t")

countsSet01 <- exprTable01[,c(3:8)]
countsSet02 <- exprTable02[,c(3:8)]
rownames(countsSet01) <- exprTable01[,c("target_id")]
rownames(countsSet02) <- exprTable02[,c("target_id")]
minCount <- 6000

# read counts
write.xlsx(exprTable01, "AvsC_DESeq2_min6000_FDR0.1_20170315.xlsx", sheetName="set01_readcounts", row.names = FALSE) 
write.xlsx(exprTable02, "AvsC_DESeq2_min6000_FDR0.1_20170315.xlsx", sheetName="set02_readcounts", row.names = FALSE, append=TRUE) 

# set01
countsAboveThr <- subset(countsSet01, rowSums(countsSet01[,c(1:6)]) > minCount)
samples <- data.frame(row.names=c("A1", "A2", "A3", "C1", "C2", "C3"), condition=as.factor(c(rep("A",3),rep("C",3))))
bckCDS2 <- DESeqDataSetFromMatrix(countData = countsAboveThr, colData=samples, design=~condition)
bckCDS2_1 <- DESeq(bckCDS2)
TCresults <- results(bckCDS2_1)
TCresults_filtered <- TCresults[ which (TCresults$padj < 0.1), ]
TCresults_filtered_table <- as.matrix(TCresults_filtered)
TCresults_filtered_table <- cbind(Row.Names = rownames(TCresults_filtered_table), TCresults_filtered_table)
rownames(TCresults_filtered_table) <- NULL
colnames(TCresults_filtered_table)[1] <- "target_id"
TCresults_filtered_table <- TCresults_filtered_table[,c(1,3,7)]
write.table(TCresults_filtered_table, file="set01_AvsC_DESeq2_min6000_FDR0.1_20170315.tsv", sep="\t", quote = FALSE, row.names = FALSE)
write.xlsx(TCresults_filtered_table, "AvsC_DESeq2_min6000_FDR0.1_20170315.xlsx", sheetName="set01_AvsC", row.names = FALSE, append=TRUE) 
# plot dendrogram
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
rld <- rlog(bckCDS2_1, blind=FALSE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
pdf("set01_AvsC_distanceheatmap_20170315.pdf")
heatmap.2(sampleDistMatrix, col = rev(colors), main = "Set 01\nSample Distances", density.info="none", trace="none", dendrogram="row")
dev.off()

# set02
countsAboveThr <- subset(countsSet02, rowSums(countsSet01[,c(1:6)]) > minCount)
samples <- data.frame(row.names=c("A1", "A2", "A3", "C1", "C2", "C3"), condition=as.factor(c(rep("A",3),rep("C",3))))
bckCDS2 <- DESeqDataSetFromMatrix(countData = countsAboveThr, colData=samples, design=~condition)
bckCDS2_1 <- DESeq(bckCDS2)
TCresults <- results(bckCDS2_1)
TCresults_filtered <- TCresults[ which (TCresults$padj < 0.1), ]
TCresults_filtered_table <- as.matrix(TCresults_filtered)
TCresults_filtered_table <- cbind(Row.Names = rownames(TCresults_filtered_table), TCresults_filtered_table)
rownames(TCresults_filtered_table) <- NULL
colnames(TCresults_filtered_table)[1] <- "target_id"
TCresults_filtered_table <- TCresults_filtered_table[,c(1,3,7)]
write.table(TCresults_filtered_table, file="set02_AvsC_DESeq2_min6000_FDR0.1_20170315.tsv", sep="\t", quote = FALSE, row.names = FALSE)
write.xlsx(TCresults_filtered_table, "AvsC_DESeq2_min6000_FDR0.1_20170315.xlsx", sheetName="set02_AvsC", row.names = FALSE, append=TRUE) 
# plot dendrogram
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
rld <- rlog(bckCDS2_1, blind=FALSE)
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
pdf("set02_AvsC_distanceheatmap_20170315.pdf")
heatmap.2(sampleDistMatrix, col = rev(colors), main = "Set 02\nSample Distances", density.info="none", trace="none", dendrogram="row")
dev.off()

# 2017-03-16
# PCA with rlog=blind
# source https://rdrr.io/bioc/DESeq2/man/rlog.html
# set #1
dds <- DESeqDataSetFromMatrix(countData = countsSet01, colData=samples, design=~condition)
rld <- rlog(dds)
dists <- dist(t(assay(rld)))
pdf("set01_complete_AvsC_clusterdendrogram_20170316.pdf")
plot(hclust(dists), main="Set #1\nTotal read counts (no filtering)")
dev.off()
# set #2
dds <- DESeqDataSetFromMatrix(countData = countsSet02, colData=samples, design=~condition)
rld <- rlog(dds)
dists <- dist(t(assay(rld)))
pdf("set02_complete_AvsC_clusterdendrogram_20170316.pdf")
plot(hclust(dists), main="Set #2\nTotal read counts (no filtering)")
dev.off()
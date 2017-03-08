setwd("/Users/barrantes/anshu/temp")
# set01
# source https://github.com/barrantesisrael/speA-2017/blob/master/dodiffexpr_parametric.R
library(DESeq2, quietly = TRUE)
exprTable <- read.csv("set01_expression_counts.tsv", header=TRUE, sep="\t")
counts <- exprTable[,c(1,3:11)]
rownames(counts) <- counts[,c("target_id")]
counts[,c("target_id")] <- NULL
minCount <- 2700
countsAboveThr <- subset(counts, rowSums(counts[,c(1:9)]) > minCount)
samples <- data.frame(row.names=c("A1", "A2", "A3", "C1", "C2", "C3", "P1", "P2", "P3"), condition=as.factor(c(rep("A",3),rep("C",3),rep("P",3))))
bckCDS2 <- DESeqDataSetFromMatrix(countData = countsAboveThr, colData=samples, design=~condition)
bckCDS2_1 <- DESeq(bckCDS2)

# Extracting transformed values
# DESeq2 manual, p.24
# "running times are shorter and more similar with blind=FALSE and 
# if the function DESeq has already been run, because then it is not
# necessary to re-estimate the dispersion values."
rld <- rlog(bckCDS2_1, blind=FALSE)

# add labels to points
# source https://support.bioconductor.org/p/79148/
library(ggplot2)
names <- colData(rld)$sample
p <- plotPCA(rld)
p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "names"), color = "black")
print(p)
# saving plot
pdf("set01.PCA.20170307.pdf")
print(p)

# set 02
exprTable <- read.csv("set02_expression_counts.tsv", header=TRUE, sep="\t")
counts <- exprTable[,c(1,3:11)]
rownames(counts) <- counts[,c("target_id")]
counts[,c("target_id")] <- NULL
minCount <- 2700
countsAboveThr <- subset(counts, rowSums(counts[,c(1:9)]) > minCount)
samples <- data.frame(row.names=c("A1", "A2", "A3", "C1", "C2", "C3", "P1", "P2", "P3"), condition=as.factor(c(rep("A",3),rep("C",3),rep("P",3))))
bckCDS2 <- DESeqDataSetFromMatrix(countData = countsAboveThr, colData=samples, design=~condition)
bckCDS2_1 <- DESeq(bckCDS2)
rld <- rlog(bckCDS2_1, blind=FALSE)
names <- colData(rld)$sample
p <- plotPCA(rld)
p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "names"), color = "black")
print(p)
pdf("set02.PCA.20170307.pdf")
print(p)
dev.off()

# modified PCAs
exprTable <- read.csv("set01_expression_counts.tsv", header=TRUE, sep="\t")
counts <- exprTable[,c(3:4,7:10)]
rownames(counts) <- exprTable[,c("target_id")]
minCount <- 2700
countsAboveThr <- subset(counts, rowSums(counts[,c(1:6)]) > minCount)
# samples <- data.frame(row.names=c("A1", "A2", "A3", "C1", "C2", "C3", "P1", "P2", "P3"), 
# condition=as.factor(c(rep("A",3),rep("C",3),rep("P",3))))
# samples <- data.frame(row.names=c("A1", "A2", "C2", "C3", "P1", "P2"), condition=as.factor(c(rep("A",2),rep("C",2),rep("P",2))))
samples <- data.frame(row.names=c("A1", "A2", "C2", "C3", "P1", "P2"), condition=as.factor(c(rep("A",2),rep("C",2),rep("P",2))))
bckCDS2 <- DESeqDataSetFromMatrix(countData = countsAboveThr, colData=samples, design=~condition)
bckCDS2_1 <- DESeq(bckCDS2)
rld <- rlog(bckCDS2_1, blind=FALSE)
names <- colData(rld)$sample
p <- plotPCA(rld)
p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "names"), color = "black")
print(p)
# pdf("set01.PCA2.20170307.pdf")
pdf("set01.PCA2.20170307.pdf")
print(p)
dev.off()

# set 2 plots
exprTable <- read.csv("set02_expression_counts.tsv", header=TRUE, sep="\t")
countsM <- exprTable[,c(4:7,9:10)]
countsN <- exprTable[,c(4:5,7:8,9:10)]
countsO <- exprTable[,c(4:5,7:8,10:11)]
rownames(countsM) <- exprTable[,c("target_id")]
rownames(countsN) <- exprTable[,c("target_id")]
rownames(countsO) <- exprTable[,c("target_id")]
# group M
countsAboveThr <- subset(countsM, rowSums(counts[,c(1:6)]) > minCount)
samples <- data.frame(row.names=c("A2", "A3", "C1", "C2", "P1", "P2"), condition=as.factor(c(rep("A",2),rep("C",2),rep("P",2))))
bckCDS2 <- DESeqDataSetFromMatrix(countData = countsAboveThr, colData=samples, design=~condition)
bckCDS2_1 <- DESeq(bckCDS2)
rld <- rlog(bckCDS2_1, blind=FALSE)
names <- colData(rld)$sample
p <- plotPCA(rld)
p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "names"), color = "black")
print(p)
pdf("set02.PCA3.20170307.pdf")
print(p)
dev.off()
# group N
countsAboveThr <- subset(countsN, rowSums(counts[,c(1:6)]) > minCount)
samples <- data.frame(row.names=c("A2", "A3", "C2", "C3", "P1", "P2"), condition=as.factor(c(rep("A",2),rep("C",2),rep("P",2))))
bckCDS2 <- DESeqDataSetFromMatrix(countData = countsAboveThr, colData=samples, design=~condition)
bckCDS2_1 <- DESeq(bckCDS2)
rld <- rlog(bckCDS2_1, blind=FALSE)
names <- colData(rld)$sample
p <- plotPCA(rld)
p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "names"), color = "black")
print(p)
pdf("set02.PCA4.20170307.pdf")
print(p)
dev.off()
# group O
countsAboveThr <- subset(countsO, rowSums(counts[,c(1:6)]) > minCount)
samples <- data.frame(row.names=c("A2", "A3", "C2", "C3", "P2", "P3"), condition=as.factor(c(rep("A",2),rep("C",2),rep("P",2))))
bckCDS2 <- DESeqDataSetFromMatrix(countData = countsAboveThr, colData=samples, design=~condition)
bckCDS2_1 <- DESeq(bckCDS2)
rld <- rlog(bckCDS2_1, blind=FALSE)
names <- colData(rld)$sample
p <- plotPCA(rld)
p <- p + geom_text(aes_string(x = "PC1", y = "PC2", label = "names"), color = "black")
print(p)
pdf("set02.PCA5.20170307.pdf")
print(p)
dev.off()

#!/usr/bin/env Rscript
# 2017-02-03 isradelacon@gmail.com
library(DESeq2, quietly = TRUE)

exprTable <- read.csv("infile", header=TRUE, sep="\t")
counts <- exprTable[,c(1,3:11)]
rownames(counts) <- counts[,c("target_id")]
counts[,c("target_id")] <- NULL
minCount <- 2700
countsAboveThr <- subset(counts, rowSums(counts[,c(1:9)]) > minCount)

samples <- data.frame(row.names=c("A1", "A2", "A3", "C1", "C2", "C3", "P1", "P2", "P3"), condition=as.factor(c(rep("A",3),rep("C",3),rep("P",3))))
bckCDS2 <- DESeqDataSetFromMatrix(countData = countsAboveThr, colData=samples, design=~condition)
bckCDS2_1 <- DESeq(bckCDS2, fitType = "mean")
ACresCtrst2 <- results(bckCDS2_1, contrast=c("condition","C","A"))
APresCtrst2 <- results(bckCDS2_1, contrast=c("condition","P","A"))
CPresCtrst2 <- results(bckCDS2_1, contrast=c("condition","P","C"))
ACresCtrst2 <- ACresCtrst2[order(ACresCtrst2$padj), ]
APresCtrst2 <- APresCtrst2[order(APresCtrst2$padj), ]
CPresCtrst2 <- CPresCtrst2[order(CPresCtrst2$padj), ]
ACresCtrst_sig2 <- ACresCtrst2[ which (ACresCtrst2$padj < 0.1), ]
APresCtrst_sig2 <- APresCtrst2[ which (APresCtrst2$padj < 0.1), ]
CPresCtrst_sig2 <- CPresCtrst2[ which (CPresCtrst2$padj < 0.1), ]
write.table(ACresCtrst_sig2, file="anshuRNAseq_AvsC_DESeq2_min2700_FDR0.1_20170130.tsv", sep="\t")
write.table(APresCtrst_sig2, file="anshuRNAseq_AvsP_DESeq2_min2700_FDR0.1_20170130.tsv", sep="\t")
write.table(CPresCtrst_sig2, file="anshuRNAseq_CvsP_DESeq2_min2700_FDR0.1_20170130.tsv", sep="\t")

sessionInfo()

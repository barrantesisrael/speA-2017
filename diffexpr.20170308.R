library(DESeq2, quietly = TRUE)

setwd("/Users/barrantes/anshu/temp")
exprTable01 <- read.csv("set01_expression_counts.tsv", header=TRUE, sep="\t")
exprTable02 <- read.csv("set02_expression_counts.tsv", header=TRUE, sep="\t")

countsSet01A <- exprTable01[,c(3:4,7:8,9:10)]
countsSet02B <- exprTable02[,c(4:5,7:8,9:10)]
rownames(countsSet01A) <- exprTable01[,c("target_id")]
rownames(countsSet02B) <- exprTable02[,c("target_id")]
minCount <- 6000

# set01, group A
countsAboveThr <- subset(countsSet01A, rowSums(countsSet01A[,c(1:6)]) > minCount)
samples <- data.frame(row.names=c("A1", "A2", "C2", "C3", "P1", "P2"), condition=as.factor(c(rep("A",2),rep("C",2),rep("P",2))))
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
write.table(ACresCtrst_sig2, file="set01_AvsC_DESeq2_min6000_FDR0.1_20170308.tsv", sep="\t")
write.table(APresCtrst_sig2, file="set01_AvsP_DESeq2_min6000_FDR0.1_20170308.tsv", sep="\t")
write.table(CPresCtrst_sig2, file="set01_CvsP_DESeq2_min6000_FDR0.1_20170308.tsv", sep="\t")

# set02, group B
countsAboveThr <- subset(countsSet02B, rowSums(countsSet02B[,c(1:6)]) > minCount)
samples <- data.frame(row.names=c("A2", "A3", "C2", "C3", "P1", "P2"), condition=as.factor(c(rep("A",2),rep("C",2),rep("P",2))))
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
write.table(ACresCtrst_sig2, file="set02_AvsC_DESeq2_min6000_FDR0.1_20170308.tsv", sep="\t")
write.table(APresCtrst_sig2, file="set02_AvsP_DESeq2_min6000_FDR0.1_20170308.tsv", sep="\t")
write.table(CPresCtrst_sig2, file="set02_CvsP_DESeq2_min6000_FDR0.1_20170308.tsv", sep="\t")

minCount <- 3000
countsAboveThr <- subset(countsSet02B, rowSums(countsSet02B[,c(1:6)]) > minCount)
samples <- data.frame(row.names=c("A2", "A3", "C2", "C3", "P1", "P2"), condition=as.factor(c(rep("A",2),rep("C",2),rep("P",2))))
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
write.table(ACresCtrst_sig2, file="set02_AvsC_DESeq2_min3000_FDR0.1_20170308.tsv", sep="\t")
write.table(APresCtrst_sig2, file="set02_AvsP_DESeq2_min3000_FDR0.1_20170308.tsv", sep="\t")
write.table(CPresCtrst_sig2, file="set02_CvsP_DESeq2_min3000_FDR0.1_20170308.tsv", sep="\t")



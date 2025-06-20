library(tidyverse)
library(DESeq2)

myCounts <- read.csv("counts.csv", row.names = 1)
colData<-read.csv("meta.csv", row.names = 1)
myCounts = myCounts[1:12]
all(colnames(myCounts) == rownames(colData))
dds <- DESeqDataSetFromMatrix(countData = myCounts, colData = colData, design = ~Group)
keep <- rowSums(counts(dds)) > 50
head(keep)
table(keep)
dds = dds[keep, ]
head(dds)
dds$Group <- relevel(dds$Group, ref = "NC")
dds <- DESeq(dds)
res <- results(dds)
summary(res)
res <- results(dds, alpha = 0.0001)
summary(res)
plotMA(res)
write.csv(as.data.frame(res), file = "resultofdeseq2.csv")


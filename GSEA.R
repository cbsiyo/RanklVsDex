library(tidyverse)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)

sigs <- read.csv("resultofdeseq2.csv", row.names = 1)
sigs = sigs[order(-sigs$stat),]
statlist <- sigs$stat
names(statlist) <- rownames(sigs)
ges <- gseGO(statlist, ont = "BP", keyType = "ENSEMBL", OrgDb = "org.Mm.eg.db", 
             minGSSize = 50, eps = 1e-10, pvalueCutoff = 1e-5,)
gesres <- as.data.frame(ges)
gesres = gesres[order(-gesres$enrichmentScore),]
write.csv(as.data.frame(ges), file = "resultofgsea.csv")
gseaplot(ges, geneSetID = "GO:0006119", title = "GESA of Oxidative Phosphorylation")
### DEseq2 analysis - R part

#Input with spike-in normalized value#

setwd("//10.10.2.7/soomin/DROSHA_N/200521_HCT116_analysis/Results")

library(DESeq2)
# library(DEGreport)
library(tidyverse)
library(pheatmap)

# Load the data normed count with spike-ins
data <- read.csv("norm_count_w_spikeIn.csv", row.names=1)
meta <- read.table("meta.txt", header=T, row.names=1)

head(data)
head(meta)

# Match the meta and data
(colnames(data) %in% rownames(meta))
(colnames(data) == rownames(meta))

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData=round(data), colData=meta, design=~sampletype)

# Sample-level QC
rld <- rlog(dds, blind=TRUE)

# PCA
plotPCA(rld, intgroup='sampletype')

# Hierarchical clustering
rld_mat <- assay(rld)
rld_cor <-cor(rld_mat)
head(rld_cor)

pheatmap(rld_cor)

# DE analysis
dds <- DESeq(dds)
sizeFactors(dds)
colSums(counts(dds))
colSums(counts(dds,normalized=TRUE))
plotDispEsts(dds)

# Define contrasts, extract results table, and shrink log2 fold changes
resultsNames(dds)
res_table_unshrunken <- results(dds, name='sampletype_sDro_vs_FL')
res_table <- lfcShrink(dds, coef='sampletype_sDro_vs_FL', type='apeglm')

head(res_table)

# MA plot with unshrunken results
plotMA(res_table_unshrunken)

# MA plot with shrunken results
plotMA(res_table)

# Save results - normed counts without spike-ins
#write.csv(as.data.frame(res_table_unshrunken %>% data.frame()), file="DE_norm_w_spikeIn_unshrink.csv")
write.csv(as.data.frame(res_table %>% data.frame()), file="220512_DE_norm_w_spikeIn_shrink.csv")

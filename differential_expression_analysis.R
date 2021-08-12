################################################################################
##  R script for differential gene expression analysis and visualization
##  Ming-an Sun, 2021-08-11
################################################################################

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(biomaRt)
require(plyr)

##### Read count table
count_table <- read.table("MN_read_count.txt", header=T, row.names = 1)
colnames(count_table) <- c("WT_1","WT_2", "KO_1","KO_2")

##### Gene annotation and filtering

gene_id <- row.names(count_table)

## Get biomaRt annotation
mart <- useMart("ensembl")
gene_inf <- getBM(
  attributes = c(
    "ensembl_gene_id","external_gene_name", "chromosome_name", "gene_biotype"
    ),  
  mart = useDataset("mmusculus_gene_ensembl",mart),
  filters ="ensembl_gene_id", values=gene_id
)

## Filter genes after add gene symbols
gene_inf_flt <- gene_inf[
  gene_inf$gene_biotype %in% c("protein_coding") & gene_inf$chromosome_name != "MT",
  ]

row.names(count_table) = mapvalues(
  row.names(count_table),
  from = gene_inf_flt$ensembl_gene_id,
  to = paste0(
    gene_inf_flt$ensembl_gene_id, ":", gene_inf_flt$external_gene_name
    )
)

count_table <- count_table[grep(":", row.names(count_table)),]

##### Differential expression analysis and visualization

## make DESeq dataset
cnt_data <- count_table[apply(count_table, 1, sum) >= 10, ]
col_data <- data.frame(Sample = c("WT", "WT", "HD", "HD"))
dds <- DESeqDataSetFromMatrix(
  countData = cnt_data, 
  colData = col_data, 
  design = ~ Sample
  )

## Normalize read counts
dds <- estimateSizeFactors(dds)
dds_nrm <- counts(dds, normalized = TRUE)
write.table(
  dds_nrm, 
  file="MN_read_count.normalized.txt", 
  sep="\t", quote=FALSE, 
  row.names = TRUE, col.names = TRUE
  )

## Correlation heatmap
dds_nrm_log10 <- log10(dds_nrm + 1)
png(file="MN_correlation_heatmap.png", width = 2400, height = 2000, res = 600)
pheatmap(
  cor(dds_nrm_log10, method = "p"), 
  display_numbers = TRUE, number_format = "%.3f"
  )
dev.off()

## PCA analysis and visualization
rld <- rlog(dds, blind = TRUE)
pca_data <- plotPCA(
  rld, intgroup = c("Sample"), 
  returnData = TRUE, ntop = 1000
  )
pca_var  <- round(100 * attr(pca_data, "percentVar"))

png(file="MN_PCA_plot.png", width=2400, height=1600, res=600)
ggplot(pca_data, aes(PC1, PC2, color = Sample)) +
  theme_bw() + xlim(-6, 6) + ylim(-6, 6) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percent.rld.var[1], "% variance")) +
  ylab(paste0("PC2: ", percent.rld.var[2], "% variance"))
dev.off()


## Call differentially expressed genes

design(dds) <- ~Sample

dds <- DESeq(dds)
res <- results(dds, contrast = c("Sample", "HD", "WT"))

up.list <- row.names(
  res[
    !is.na(res$log2FoldChange) & 
      !is.na(res$padj) & 
      res$log2FoldChange > 1 
    & res$padj < 0.05,
    ]
  )

down.list <- row.names(
  res[
    !is.na(res$log2FoldChange) & 
      !is.na(res$padj) & 
      res$log2FoldChange < -1 & 
      res$padj < 0.05,
    ]
  )

## Exclude imprinted genes which usually have heterogenous expression among colonies
imprinted_genes <- read.table("Mouse.geneImprint.full.txt")[,1]

up.list   <-   up.list[!sub("^.*:", "",   up.list) %in% imprinted_genes]
down.list <- down.list[!sub("^.*:", "", down.list) %in% imprinted_genes]

res.flt <- res[!sub("^.*:", "", row.names(res)) %in% imprinted_genes,]
res.srt <- res.flt[order(res.flt$padj),]

de.list <- unique(c(up.list, down.list))

## Save DE results as table
write.table(
  res.srt, 
  file = "MN.DEtable.txt", 
  row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE
)

## Generate MA-plot for DEGs
png(file="MN_MA_plot.png", width = 3000, height = 2600, res = 600)
par(mar=c(4,4,1,1), mfrow = c(1,1))
plot(
  log10(res$baseMean+1), res$log2FoldChange, 
  xlab = "Average read counts", ylab = "log2(HD/WT)", 
  xlim = c(0,5), ylim = c(-5,5),
  xaxt = "n", cex = .2, pch = 19, col = alpha('black',0.5)
  )
axis(1, at = 0:5, labels = 10^(0:5))
abline(h = 0, lwd = 3, col = alpha('red', 0.5))
points(
  log10(res$baseMean[row.names(res) %in% up.list] + 1), 
  res$log2FoldChange[row.names(res) %in% up.list], 
  cex = 0.5, pch = 19, col = 'red'
  )
points(
  log10(res$baseMean[row.names(res) %in% down.list] + 1), 
  res$log2FoldChange[row.names(res) %in% down.list], 
  cex = 0.5, pch = 19, col = 'green'
  )
dev.off()

# Generate heatmap for DEGs
de.cnt.nrm <- dds.cnt.nrm[row.names(dds.cnt.nrm) %in% de.list, ]
de.cnt.ann <- data.frame(
  row.names = de.list,
  Expression = factor(
    rep("Unchanged", length(de.list)),
    levels = c("Up-regulated", "Down-regulated", "Unchanged")
    )
)

de.cnt.ann$Expression[row.names(de.cnt.ann) %in%   up.list] = "Up-regulated"
de.cnt.ann$Expression[row.names(de.cnt.ann) %in% down.list] = "Down-regulated"
de.cnt.ann$Expression = factor(
  de.cnt.ann$Expression, 
  levels = c("Up-regulated", "Down-regulated", "Unchanged")
  )

png(file = "Fig.HD2WT_DE_heatmap.png", width = 3200, height = 3000, res = 600)
pheatmap(
  log10(de.cnt.nrm+1), scale = "row", 
  show_rownames = FALSE, annotation_row = de.cnt.ann,
  border_color = NA,
  annotation_colors = list(
    Expression = c(
      "Up-regulated" = "Red", 
      "Down-regulated" = "Green", 
      "Unchanged" = "Grey"
      )
  )
)
dev.off()

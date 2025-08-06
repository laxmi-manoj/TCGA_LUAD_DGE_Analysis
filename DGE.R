# Lung Cancer RNA-seq DEG Analysis with DESeq2

# =========================
# CONFIGURATION
# =========================
log2fc_thresh    <- 1
padj_thresh      <- 0.05
top_n_label      <- 5
top_n_heatmap    <- 30
top_n_bar        <- 10
dirs             <- c("data", "results", "figures")
pkgs             <- c(
  "BiocManager","DESeq2","ggplot2","pheatmap","clusterProfiler",
  "org.Hs.eg.db","EnhancedVolcano","RColorBrewer","AnnotationDbi",
  "tidyverse","GEOquery","readr","apeglm","forcats","TCGAbiolinks",
  "scales","biomaRt"
)

# =========================
# UTILITIES
# =========================
# install & load
for (pkg in pkgs) {
  if (!require(pkg, character.only=TRUE)) {
    if (pkg=="BiocManager") install.packages("BiocManager")
    else BiocManager::install(pkg, ask=FALSE)
    library(pkg, character.only=TRUE)
  }
}
library(dplyr)

# create folders
lapply(dirs, dir.create, recursive=TRUE, showWarnings=FALSE)

# helper to save plots
save_plot <- function(plot, file) {
  ggsave(file, plot=plot, width=10, height=8, dpi=300)
}

# =========================
# 1. DATA DOWNLOAD & PREP
# =========================
query  <- GDCquery(
  project="TCGA-LUAD", data.category="Transcriptome Profiling",
  data.type="Gene Expression Quantification",
  workflow.type="STAR - Counts",
  sample.type=c("Primary Tumor","Solid Tissue Normal")
)
GDCdownload(query, method="api", files.per.chunk=20)
data   <- GDCprepare(query)
counts <- assay(data)
coldata<- colData(data)
metadata <- data.frame(condition=ifelse(coldata$sample_type=="Primary Tumor","Tumor","Normal"))
rownames(metadata) <- colnames(counts)
write.csv(counts, "data/counts.csv", row.names=TRUE)
write.csv(metadata, "data/metadata.csv", row.names=TRUE)

# =========================
# 2. DESEQ2 ANALYSIS
# =========================
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~condition)
dds <- dds[rowSums(counts(dds))>10, ]
dds <- DESeq(dds)

res <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")

res_df <- as.data.frame(res) %>%
  mutate(
    gene_id       = rownames(.),
    ensembl_clean = sub("\\.\\d+$","",gene_id),
    gene_name     = mapIds(org.Hs.eg.db, keys=ensembl_clean, column="SYMBOL", keytype="ENSEMBL"),
    entrez_id     = mapIds(org.Hs.eg.db, keys=ensembl_clean, column="ENTREZID", keytype="ENSEMBL"),
    Significance  = if_else(padj<padj_thresh & abs(log2FoldChange)>log2fc_thresh, "Significant", "Not Significant")
  ) %>%
  filter(!is.na(log2FoldChange), !is.na(padj))

# save full results
write.csv(res_df, "results/all_deseq2_results.csv", row.names=FALSE)

# =========================
# 3. SIGNIFICANT DEGS
# =========================
sig <- res_df %>%
  filter(padj<padj_thresh, abs(log2FoldChange)>log2fc_thresh) %>%
  arrange(desc(abs(log2FoldChange)))

sig_out <- dplyr::select(sig,
                         gene_id,
                         ensembl_clean,
                         entrez_id,
                         biomart_id = ensembl_clean,
                         gene_name,
                         log2FoldChange,
                         pvalue,
                         padj
)

write.csv(sig_out, "results/significant_degs.csv", row.names=FALSE)

# =========================
# 4. QC & PLOTS
# =========================
vsd <- vst(dds, blind=FALSE)

## PCA
pcaData    <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p_pca <- ggplot(pcaData, aes(PC1,PC2,color=condition)) +
  geom_point(size=5,alpha=0.8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw(base_size=16) + ggtitle("PCA Plot")
save_plot(p_pca, "figures/pca_plot.png")

## Heatmap
top30   <- sig %>% slice_min(padj, n=top_n_heatmap)
mat     <- assay(vsd)[top30$gene_id, ]
rownames(mat) <- coalesce(top30$gene_name, top30$gene_id)
pheatmap(mat,
         annotation_col=data.frame(Condition=metadata$condition, row.names=rownames(metadata)),
         show_colnames=FALSE, scale="row", fontsize_row=10,
         filename="figures/heatmap_top30.png", width=10, height=8)

## Volcano
res_df$label <- coalesce(res_df$gene_name, res_df$gene_id)
to_label     <- res_df %>% arrange(padj) %>% slice_head(n=top_n_label) %>% pull(label)
volcano_plot <- EnhancedVolcano(
  res_df, lab=res_df$label, selectLab=to_label,
  x='log2FoldChange', y='padj',
  pCutoff=padj_thresh, FCcutoff=log2fc_thresh,
  title='Volcano Plot: Tumor vs Normal'
)
save_plot(volcano_plot, "figures/volcano_plot.png")

## Barplot
top_up   <- slice_max(sig, log2FoldChange, n=top_n_bar)
top_down <- slice_min(sig, log2FoldChange, n=top_n_bar)
top_genes<- bind_rows(top_up, top_down) %>%
  mutate(gene=coalesce(gene_name, gene_id))
p_bar <- ggplot(top_genes, aes(x=fct_reorder(gene,log2FoldChange), y=log2FoldChange, fill=log2FoldChange>0)) +
  geom_col() + coord_flip() + theme_bw(base_size=16) +
  xlab("Gene") + ylab("Log2 Fold Change") + ggtitle("Top 10 Up/Down Genes")
save_plot(p_bar, "figures/barplot_top_genes.png")

# =========================
# 5. FUNCTIONAL ENRICHMENT
# =========================
entrez_list <- na.omit(unique(sig$entrez_id))

## GO
ego        <- enrichGO(
  gene=entrez_list, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
  ont="BP", pAdjustMethod="BH", pvalueCutoff=1, qvalueCutoff=1, readable=TRUE
)
go_results <- as.data.frame(ego)
write.csv(go_results, "results/go_enrichment_results.csv", row.names=FALSE)
barplot(ego, showCategory=10); save_plot(last_plot(), "figures/go_enrichment.png")

## KEGG
ekk         <- enrichKEGG(gene=entrez_list, organism='hsa', pvalueCutoff=1, qvalueCutoff=1)
kegg_results<- as.data.frame(ekk)
write.csv(kegg_results, "results/kegg_enrichment_results.csv", row.names=FALSE)
barplot(ekk, showCategory=10); save_plot(last_plot(), "figures/kegg_enrichment.png")

# =========================
# 6. SUMMARY & REPRODUCTION
# =========================
writeLines(capture.output(sessionInfo()), "results/session_info.txt")

cat(paste0(
  "\nSummary:\n",
  "Tumor samples: ", sum(metadata$condition=="Tumor"), "\n",
  "Normal samples: ", sum(metadata$condition=="Normal"), "\n",
  "Significant DEGs: ", nrow(sig), "\n",
  "GO terms enriched: ", nrow(go_results), "\n",
  "KEGG pathways enriched: ", nrow(kegg_results), "\n"
))

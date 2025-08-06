# TCGA_LUAD_DGE

# Differential Gene Expression Analysis in Lung Adenocarcinoma (LUAD)

Project by Abel George and Laxmi CM

This repository contains an end-to-end RNA-seq analysis pipeline for identifying differentially expressed genes (DEGs) in **lung adenocarcinoma (LUAD)** using data from **The Cancer Genome Atlas (TCGA)**. The analysis was implemented in R and includes downstream functional enrichment and protein–protein interaction (PPI) network analysis using STRING and Cytoscape.

## 📁 Project Structure
├── DGE.R # Main R script for the complete workflow
├── results/ # Output files (DEG tables, enrichment results, etc.)
├── figures/ # Plots (PCA, heatmap, volcano, barplots)
├── network/ # STRING/Cytoscape node and edge tables, network export
└── README.md # Project documentation (this file)

---

## 🧪 Analysis Overview

1. **Data Source**  
   - TCGA-LUAD RNA-seq count data (via `TCGAbiolinks`)
   - Primary Tumor (n = 539) vs. Normal Tissue (n = 59)

2. **Tools & Packages**  
   - `DESeq2`, `apeglm`, `ggplot2`, `EnhancedVolcano`, `pheatmap`,  
     `clusterProfiler`, `biomaRt`, `org.Hs.eg.db`, `STRINGdb`, `Cytoscape`

3. **Pipeline Steps**
   - Download & preprocess TCGA data
   - Normalize and filter counts
   - Identify DEGs (|log2FC| > 1, FDR < 0.05)
   - Visualize DEGs (PCA, heatmap, volcano, barplot)
   - Functional enrichment (GO, KEGG)
   - PPI network construction via STRING (confidence > 0.7)  
   - Cytoscape-based visualization (nodes colored by log2FC, edges by confidence)

---

## 📊 Key Results

- **DEGs Identified**:  
  14,789 genes (11,712 upregulated, 3,077 downregulated)

- **Top DEGs by significance**:  
  FAM83A, PYCR1, OTUD1, AFAP1-AS1, EPAS1

- **Top DEGs by log2FC**:  
  PSG1, TFF2, MAGEA9, REG4, DEFA5 (up)  
  SLC6A4, AGER, FABP4, SFTPC, CSF3 (down)

- **Top GO Terms**:  
  - Extracellular matrix organization  
  - Humoral immune response  
  - GPCR signaling  
  - Chromatin localization

- **Top KEGG Pathways**:  
  - Neuroactive ligand-receptor interaction  
  - Systemic lupus erythematosus  
  - Neutrophil extracellular trap formation  
  - Hormone signaling

- **PPI Network**:  
  - Top 200 upregulated & 200 downregulated genes submitted to STRING  
  - Confidence cutoff = 0.7  
  - Network visualized in Cytoscape with styling based on log2FC and edge confidence

---

## 📂 Outputs

- **CSV tables**:  
  - `results/DEG_upregulated.csv`, `DEG_downregulated.csv`  
  - `results/GO_enrichment.csv`, `KEGG_enrichment.csv`

- **Plots** (in `figures/`):  
  - `volcano_plot.png`  
  - `pca_plot.png`  
  - `heatmap_top30.png`  
  - `barplot_top_genes.png`  
  - `go_barplot.png`, `kegg_barplot.png`  

- **Network Data** (in `network/`):  
  - `STRING network default node.csv`  
  - `STRING network default edge.csv`
  - `ppi_network_image.png`
    
---

## ▶️ How to Run

```r
source("DGE.R")



# scRNA-spermatogenesis
This repository contains R scripts for the analysis of mouse single-cell RNA sequencing (scRNA-seq) data from germ cells to investigate transcriptomic changes throughout in vivo spermatogenesis.

## Analysis
The repository includes R scripts for the following analyses:

1. **Data acquisition and preprocessing**
   - Quality control (QC) of each scRNA-seq dataset.
   - Identification and extraction of germ cells from each dataset.

2. **Data integration**
   - Integration of all datasets into a single Seurat object containing only germ cells.

3. **Dimensionality reduction and clustering**
   - UMAP visualization.
   - Unsupervised clustering of germ cells.

4. **Cell type annotation**
   - Identification of major germ cell populations:
     - Spermatogonia
     - Spermatocytes
     - Round spermatids
     - Elongated spermatids

5. **Gene co-expression analysis**
   - Construction of gene co-expression networks using hdWGCNA.

6. **Lipid metabolism analysis**
   - Co-expression analysis of lipid-related genes using hdWGCNA.
   - Identification of differentially expressed lipid-related genes across germ cell populations.

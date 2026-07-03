#Cargar librerías
library (scDblFinder)
library (clusterProfiler)
library (Seurat)
library (Matrix)
library (hdf5r)
library (ggplot2)
library (tidyverse)
library (dplyr)
library(patchwork)
library (harmony)

############GSM3744441 - E18.5##########
GSM3744441_E18.5_barcodes <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM3744441_E18.5_Rep1/GSM3744441_E18.5_Rep1_barcodes.tsv.gz"
GSM3744441_E18.5_feature <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM3744441_E18.5_Rep1/GSM3744441_E18.5_Rep1_features.tsv.gz"
GSM3744441_E18.5_matrix <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM3744441_E18.5_Rep1/GSM3744441_E18.5_Rep1_matrix.mtx.gz"
#1- Creación del objeto Seurat
counts_GSM3744441 <- readMM(GSM3744441_E18.5_matrix)
barcodes_GSM3744441 <- read.table(GSM3744441_E18.5_barcodes, stringsAsFactors=F)[,1]
features_GSM3744441 <- read.csv(GSM3744441_E18.5_feature, stringsAsFactors=F, sep="\t", header=F)
rownames(counts_GSM3744441) <- make.unique(features_GSM3744441[,2])
colnames(counts_GSM3744441) <- barcodes_GSM3744441
GSM3744441 <- CreateSeuratObject(counts_GSM3744441, 
                                 min.cells = 3, 
                                 min.features = 1000, 
                                 project="E18.5")
rm (features_GSM3744441, barcodes_GSM3744441, counts_GSM3744441)
rm (GSM3744441_E18.5_barcodes, GSM3744441_E18.5_feature, GSM3744441_E18.5_matrix)
dim (GSM3744441)  
head (GSM3744441@meta.data)
#2- Análisis de calidad
GSM3744441[["percent.mt"]] <- PercentageFeatureSet(GSM3744441, pattern = "^mt[-\\.]")
head (GSM3744441@meta.data)
plot1 <- FeatureScatter(GSM3744441, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(GSM3744441, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2
GSM3744441 <- subset(GSM3744441, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
plot3 <- FeatureScatter(GSM3744441, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(GSM3744441, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
rm (plot1, plot2, plot3, plot4)
dim (GSM3744441)  
#3- Eliminar multipletes (droplets)
#Creamos un archivo tipo 'Single Cell'
GSM3744441_sce <- as.SingleCellExperiment(GSM3744441)
#Corremos el scrip para la busqueda de dobletes y visualizamos el resultado
GSM3744441_sce <- scDblFinder(GSM3744441_sce, samples = "orig.ident") #Tarda unos minutos. 
head (GSM3744441_sce)
table(GSM3744441_sce$scDblFinder.class)
#Exploramos el resultado y lo incorporamos al objeto Seurat:
meta_scdblfinder_GSM3744441 <- GSM3744441_sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
head(meta_scdblfinder_GSM3744441)
rownames(meta_scdblfinder_GSM3744441) <- GSM3744441_sce@colData@rownames
GSM3744441 <- AddMetaData(object = GSM3744441, metadata = meta_scdblfinder_GSM3744441 %>% dplyr::select('scDblFinder.class'))
VlnPlot(GSM3744441, split.by = "scDblFinder.class", #ver los de sampleID
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
GSM3744441 <- subset(GSM3744441, subset = scDblFinder.class == 'singlet')
rm (meta_scdblfinder_GSM3744441, GSM3744441_sce)
dim (GSM3744441)
GSM3744441 <- subset(GSM3744441, subset = nCount_RNA <= 30000)
#4- Filtrar celulas germinales
GSM3744441 <- NormalizeData(GSM3744441)
GSM3744441 <- FindVariableFeatures(GSM3744441, nfeatures = 2000)
GSM3744441 <- ScaleData(GSM3744441)
GSM3744441 <-  RunPCA(GSM3744441)
ElbowPlot(GSM3744441)
GSM3744441 <-  RunUMAP(GSM3744441, dims = 1:6) #pUEDO SEPARAR MEJOR LAS GERMINALES
UMAPPlot(GSM3744441)
FeaturePlot(GSM3744441, features= "Ddx4") #Germinal
FeaturePlot(GSM3744441, features= "Dazl") #Germinal
FeaturePlot(GSM3744441, features= "Sox9") #Sertoli
FeaturePlot(GSM3744441, features= "Amh") #Sertoli
FeaturePlot(GSM3744441, features= "Esam") #Endotelial
FeaturePlot(GSM3744441, features= "Pecam1") #Endotelial
FeaturePlot(GSM3744441, features= "Lyz2") #Macrofago
FeaturePlot(GSM3744441, features= "Cd68") #Macrofago
FeaturePlot(GSM3744441, features= "Cnmd") #Inmune
FeaturePlot(GSM3744441, features= "Cyp11a1") #lEY
FeaturePlot(GSM3744441, features= "Igf1") #Estromal
FeaturePlot(GSM3744441, features= "Acta2") #peritubular
GSM3744441 <-  FindNeighbors(GSM3744441, dims = 1:6)
GSM3744441 <- FindClusters(GSM3744441, resolution = c(0.2)) #11 COMUNIDADES
DimPlot(GSM3744441, reduction = "umap", label = TRUE) #RESOLUTION.0.2 == CLUSTER 10
VlnPlot(GSM3744441, features = "Ddx4")
VlnPlot(GSM3744441, features = "Mcm2")
VlnPlot(GSM3744441, features = "Dazl") #Verifico
GSM3744441_germinales <- subset (GSM3744441, ident= 9)
dim (GSM3744441_germinales) 

#GSM3744441 <- SCTransform(GSM3744441, verbose = FALSE)

############GSM3744442 - E18.5##########
GSM3744442_E18.5_barcodes <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM3744442_E18.5_Rep2/GSM3744442_E18.5_Rep2_barcodes.tsv.gz"
GSM3744442_E18.5_feature <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM3744442_E18.5_Rep2/GSM3744442_E18.5_Rep2_features.tsv.gz"
GSM3744442_E18.5_matrix <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM3744442_E18.5_Rep2/GSM3744442_E18.5_Rep2_matrix.mtx.gz"
#1- Creación del objeto Seurat
counts_GSM3744442 <- readMM(GSM3744442_E18.5_matrix)
barcodes_GSM3744442 <- read.table(GSM3744442_E18.5_barcodes, stringsAsFactors=F)[,1]
features_GSM3744442 <- read.csv(GSM3744442_E18.5_feature, stringsAsFactors=F, sep="\t", header=F)
rownames(counts_GSM3744442) <- make.unique(features_GSM3744442[,2])
colnames(counts_GSM3744442) <- barcodes_GSM3744442
GSM3744442 <- CreateSeuratObject(counts_GSM3744442, min.cells = 3, min.features = 1000, project="GSM3744442_E18.5")
rm (features_GSM3744442, barcodes_GSM3744442, counts_GSM3744442)
rm (GSM3744442_E18.5_barcodes, GSM3744442_E18.5_feature, GSM3744442_E18.5_matrix)
dim (GSM3744442)  
head (GSM3744442@meta.data)
#2- Análisis de calidad
GSM3744442[["percent.mt"]] <- PercentageFeatureSet(GSM3744442, pattern = "^mt[-\\.]")
head (GSM3744442@meta.data)
plot1 <- FeatureScatter(GSM3744442, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(GSM3744442, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2
GSM3744442 <- subset(GSM3744442, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
plot3 <- FeatureScatter(GSM3744442, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(GSM3744442, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
rm (plot1, plot2, plot3, plot4)
dim (GSM3744442)  
#3- Eliminar multipletes (droplets)
#Creamos un archivo tipo 'Single Cell'
GSM3744442_sce <- as.SingleCellExperiment(GSM3744442)
#Corremos el scrip para la busqueda de dobletes y visualizamos el resultado
GSM3744442_sce <- scDblFinder(GSM3744442_sce, samples = "orig.ident") #Tarda unos minutos. 
head (GSM3744442_sce)
table(GSM3744442_sce$scDblFinder.class)
#Exploramos el resultado y lo incorporamos al objeto Seurat:
meta_scdblfinder_GSM3744442 <- GSM3744442_sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
head(meta_scdblfinder_GSM3744442)
rownames(meta_scdblfinder_GSM3744442) <- GSM3744442_sce@colData@rownames
GSM3744442 <- AddMetaData(object = GSM3744442, metadata = meta_scdblfinder_GSM3744442 %>% dplyr::select('scDblFinder.class'))
VlnPlot(GSM3744442, split.by = "scDblFinder.class", #ver los de sampleID
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
GSM3744442 <- subset(GSM3744442, subset = scDblFinder.class == 'singlet')
rm (meta_scdblfinder_GSM3744442, GSM3744442_sce)
dim (GSM3744442)
GSM3744442 <- subset(GSM3744442, subset = nCount_RNA <= 30000)
#4- Filtrar celulas germinales
GSM3744442 <- NormalizeData(GSM3744442)
GSM3744442 <- FindVariableFeatures(GSM3744442, nfeatures = 2000)
GSM3744442 <- ScaleData(GSM3744442)
GSM3744442 <-  RunPCA(GSM3744442)
ElbowPlot(GSM3744442)
GSM3744442 <-  RunUMAP(GSM3744442, dims = 1:6) #pUEDO SEPARAR MEJOR LAS GERMINALES
UMAPPlot(GSM3744442)
FeaturePlot(GSM3744442, features= "Ddx4") #Germinal
FeaturePlot(GSM3744442, features= "Dazl") #Germinal
FeaturePlot(GSM3744442, features= "Sox9") #Sertoli
FeaturePlot(GSM3744442, features= "Amh") #Sertoli
FeaturePlot(GSM3744442, features= "Esam") #Endotelial
FeaturePlot(GSM3744442, features= "Pecam1") #Endotelial
FeaturePlot(GSM3744442, features= "Lyz2") #Macrofago
FeaturePlot(GSM3744442, features= "Cd68") #Macrofago
FeaturePlot(GSM3744442, features= "Cnmd") #Inmune
FeaturePlot(GSM3744442, features= "Cyp11a1") #lEY
FeaturePlot(GSM3744442, features= "Igf1") #Estromal
FeaturePlot(GSM3744442, features= "Acta2") #peritubular
GSM3744442 <-  FindNeighbors(GSM3744442, dims = 1:6)
GSM3744442 <- FindClusters(GSM3744442, resolution = c(0.2)) #11 COMUNIDADES
DimPlot(GSM3744442, reduction = "umap", label = TRUE) #RESOLUTION.0.2 == CLUSTER 10
VlnPlot(GSM3744442, features = "Ddx4")
VlnPlot(GSM3744442, features = "Dazl") #Verifico
GSM3744442_germinales <- subset (GSM3744442, ident= 8)
dim (GSM3744442_germinales) #98 celulas!!!!
FeaturePlot(GSM3744442_germinales, features= "Star")
VlnPlot(GSM3744442_germinales, features= "Star")
FeaturePlot(GSM3744442_germinales, features= "nFeature_RNA")
VlnPlot(GSM3744442_germinales, features= "nFeature_RNA")

############GSM3744443 - 2 días##########
GSM3744443_P2_barcodes <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM3744443_P2_Rep1/GSM3744443_P2_Rep1_barcodes.tsv.gz"
GSM3744443_P2_feature <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM3744443_P2_Rep1/GSM3744443_P2_Rep1_features.tsv.gz"
GSM3744443_P2_matrix <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM3744443_P2_Rep1/GSM3744443_P2_Rep1_matrix.mtx.gz"
#1- Creación del objeto Seurat
counts_GSM3744443 <- readMM(GSM3744443_P2_matrix)
barcodes_GSM3744443 <- read.table(GSM3744443_P2_barcodes, stringsAsFactors=F)[,1]
features_GSM3744443 <- read.csv(GSM3744443_P2_feature, stringsAsFactors=F, sep="\t", header=F)
rownames(counts_GSM3744443) <- make.unique(features_GSM3744443[,2])
colnames(counts_GSM3744443) <- barcodes_GSM3744443
GSM3744443 <- CreateSeuratObject(counts_GSM3744443, min.cells = 3, min.features = 1000, project="GSM3744443_2days")
rm (features_GSM3744443, barcodes_GSM3744443, counts_GSM3744443)
rm (GSM3744443_P2_barcodes, GSM3744443_P2_feature, GSM3744443_P2_matrix)

dim (GSM3744443)  
head (GSM3744443@meta.data)

#2- Análisis de calidad
GSM3744443[["percent.mt"]] <- PercentageFeatureSet(GSM3744443, pattern = "^mt[-\\.]")
head (GSM3744443@meta.data)

plot1 <- FeatureScatter(GSM3744443, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(GSM3744443, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2

GSM3744443 <- subset(GSM3744443, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
plot3 <- FeatureScatter(GSM3744443, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(GSM3744443, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
rm (plot1, plot2, plot3, plot4)

dim (GSM3744443)  

#3- Eliminar multipletes (droplets)
#Creamos un archivo tipo 'Single Cell'
GSM3744443_sce <- as.SingleCellExperiment(GSM3744443)
#Corremos el scrip para la busqueda de dobletes y visualizamos el resultado
GSM3744443_sce <- scDblFinder(GSM3744443_sce, samples = "orig.ident") #Tarda unos minutos. 
head (GSM3744443_sce)
table(GSM3744443_sce$scDblFinder.class)
#Exploramos el resultado y lo incorporamos al objeto Seurat:
meta_scdblfinder_GSM3744443 <- GSM3744443_sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
head(meta_scdblfinder_GSM3744443)
rownames(meta_scdblfinder_GSM3744443) <- GSM3744443_sce@colData@rownames
GSM3744443 <- AddMetaData(object = GSM3744443, metadata = meta_scdblfinder_GSM3744443 %>% dplyr::select('scDblFinder.class'))
VlnPlot(GSM3744443, split.by = "scDblFinder.class", #ver los de sampleID
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
GSM3744443 <- subset(GSM3744443, subset = scDblFinder.class == 'singlet')
rm (meta_scdblfinder_GSM3744443, GSM3744443_sce)
dim (GSM3744443)
GSM3744443 <- subset(GSM3744443, subset = nCount_RNA < 30000)
df_GSM3744443  <- data.frame (GSM3744443@meta.data)
df_GSM3744443$QC_pass <- (df_GSM3744443$nCount_RNA < 30000 & df_GSM3744443$nFeature_RNA < 7000 & df_GSM3744443$percent.mt < 10 & df_GSM3744443$scDblFinder.class == "singlet")


GSM3744443 <- NormalizeData(GSM3744443)
GSM3744443 <- FindVariableFeatures(GSM3744443, nfeatures = 2000)
GSM3744443 <- ScaleData(GSM3744443)
GSM3744443 <-  RunPCA(GSM3744443)
ElbowPlot(GSM3744443)

GSM3744443 <-  RunUMAP(GSM3744443, dims = 1:6) #pUEDO SEPARAR MEJOR LAS GERMINALES
UMAPPlot(GSM3744443)
FeaturePlot(GSM3744443, features= "Ddx4") #Germinal
FeaturePlot(GSM3744443, features= "Dazl") #Germinal
FeaturePlot(GSM3744443, features= "Sox9") #Sertoli
FeaturePlot(GSM3744443, features= "Amh") #Sertoli
FeaturePlot(GSM3744443, features= "Esam") #Endotelial
FeaturePlot(GSM3744443, features= "Pecam1") #Endotelial
FeaturePlot(GSM3744443, features= "Lyz2") #Macrofago
FeaturePlot(GSM3744443, features= "Cd68") #Macrofago
FeaturePlot(GSM3744443, features= "Cnmd") #Inmune
FeaturePlot(GSM3744443, features= "Cyp11a1") #lEY
FeaturePlot(GSM3744443, features= "Igf1") #Estromal
FeaturePlot(GSM3744443, features= "Acta2") #peritubular
GSM3744443 <-  FindNeighbors(GSM3744443, dims = 1:6)
GSM3744443 <- FindClusters(GSM3744443, resolution = c(0.3)) #11 COMUNIDADES
DimPlot(GSM3744443, reduction = "umap", label = TRUE) #RESOLUTION.0.2 == CLUSTER 10
VlnPlot(GSM3744443, features = "Ddx4")
VlnPlot(GSM3744443, features = "Dazl")
#cl_markers_GSM3744443 <- FindAllMarkers(GSM3744443,
 #                                       features = VariableFeatures(GSM3744443))
#cl_markers_GSM3744443 <- cl_markers_GSM3744443[cl_markers_GSM3744443$avg_log2FC >= 1.5 & cl_markers_GSM3744443$p_val_adj <= 0.05,]
#df_GSM3744443 <- GSM3744443@meta.data
#df = data.frame(RNA_snn_res.0.3 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
 #               cell_type = c("Sertoli cells", "Stromal cells", "Sertoli cells", "Sertoli cells", "Stromal cells", "Myoid cells", "Stromal cells", "Myoid cells",  "Germ cells", "Leydig", "Endothelial cells", "Endothelial cells"))
#df_GSM3744443$cells <- rownames(df_GSM3744443)
#df_GSM3744443 <- data.frame (merge (df_GSM3744443, df, by= "RNA_snn_res.0.3"))
#GSM3744443_UMAP <- data.frame (Embeddings(GSM3744443, reduction = "umap"))
#GSM3744443_UMAP$cells <- rownames(GSM3744443_UMAP)

#df_GSM3744443 <- merge(df_GSM3744443, GSM3744443_UMAP, by= "cells")
#write.table(df_GSM3744443, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM3744443")
#write.table(cl_markers_GSM3744443, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM3744443_markers")

GSM3744443_germinales <- subset (GSM3744443, ident= 8)
dim (GSM3744443_germinales)

############GSM374444 - 2 días##########
GSM3744444_P2_barcodes <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM3744444_P2_Rep2/GSM3744444_P2_Rep2_barcodes.tsv.gz"
GSM3744444_P2_feature <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM3744444_P2_Rep2/GSM3744444_P2_Rep2_features.tsv.gz"
GSM3744444_P2_matrix <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM3744444_P2_Rep2/GSM3744444_P2_Rep2_matrix.mtx.gz"

#1- Creación del objeto Seurat
counts_GSM3744444 <- readMM(GSM3744444_P2_matrix)
barcodes_GSM3744444 <- read.table(GSM3744444_P2_barcodes, stringsAsFactors=F)[,1]
features_GSM3744444 <- read.csv(GSM3744444_P2_feature, stringsAsFactors=F, sep="\t", header=F)
rownames(counts_GSM3744444) <- make.unique(features_GSM3744444[,2])
colnames(counts_GSM3744444) <- barcodes_GSM3744444
GSM3744444 <- CreateSeuratObject(counts_GSM3744444, min.cells = 3, min.features = 1000, project="GSM3744444_2days")
rm (features_GSM3744444, barcodes_GSM3744444, counts_GSM3744444)
rm (GSM3744444_P2_barcodes, GSM3744444_P2_feature, GSM3744444_P2_matrix)

dim (GSM3744444)  
head (GSM3744444@meta.data)

#2- Análisis de calidad
GSM3744444[["percent.mt"]] <- PercentageFeatureSet(GSM3744444, pattern = "^mt[-\\.]")
head (GSM3744444@meta.data)

plot1 <- FeatureScatter(GSM3744444, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(GSM3744444, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2

GSM3744444 <- subset(GSM3744444, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
plot3 <- FeatureScatter(GSM3744444, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(GSM3744444, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
rm (plot1, plot2, plot3, plot4)

dim (GSM3744444)  

#3- Eliminar multipletes (droplets)
#Creamos un archivo tipo 'Single Cell'
GSM3744444_sce <- as.SingleCellExperiment(GSM3744444)
#Corremos el scrip para la busqueda de dobletes y visualizamos el resultado
GSM3744444_sce <- scDblFinder(GSM3744444_sce, samples = "orig.ident") #Tarda unos minutos. 
head (GSM3744444_sce)
table(GSM3744444_sce$scDblFinder.class)
#Exploramos el resultado y lo incorporamos al objeto Seurat:
meta_scdblfinder_GSM3744444 <- GSM3744444_sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
head(meta_scdblfinder_GSM3744444)
rownames(meta_scdblfinder_GSM3744444) <- GSM3744444_sce@colData@rownames
GSM3744444 <- AddMetaData(object = GSM3744444, metadata = meta_scdblfinder_GSM3744444 %>% dplyr::select('scDblFinder.class'))
VlnPlot(GSM3744444, split.by = "scDblFinder.class", #ver los de sampleID
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
GSM3744444 <- subset(GSM3744444, subset = scDblFinder.class == 'singlet')
rm (meta_scdblfinder_GSM3744444, GSM3744444_sce)
dim (GSM3744444)
GSM3744444 <- subset(GSM3744444, subset = nCount_RNA < 30000)
df_GSM3744444  <- data.frame (GSM3744444@meta.data)
df_GSM3744444$QC_pass <- (df_GSM3744444$nCount_RNA < 30000 & df_GSM3744444$nFeature_RNA < 7000 & df_GSM3744444$percent.mt < 10 & df_GSM3744444$scDblFinder.class == "singlet")
write.table(df_GSM3744444, file= "~/Figuras_Paper 2026/Planillas QC/GSM3744444")

GSM3744444 <- NormalizeData(GSM3744444)
GSM3744444 <- FindVariableFeatures(GSM3744444, nfeatures = 2000)
GSM3744444 <- ScaleData(GSM3744444)
GSM3744444 <-  RunPCA(GSM3744444)
ElbowPlot(GSM3744444)

GSM3744444 <-  RunUMAP(GSM3744444, dims = 1:6) #pUEDO SEPARAR MEJOR LAS GERMINALES
UMAPPlot(GSM3744444)
FeaturePlot(GSM3744444, features= "Ddx4") #Germinal
FeaturePlot(GSM3744444, features= "Dazl") #Germinal
FeaturePlot(GSM3744444, features= "Sox9") #Sertoli
FeaturePlot(GSM3744444, features= "Amh") #Sertoli
FeaturePlot(GSM3744444, features= "Esam") #Endotelial
FeaturePlot(GSM3744444, features= "Pecam1") #Endotelial
FeaturePlot(GSM3744444, features= "Lyz2") #Macrofago
FeaturePlot(GSM3744444, features= "Cd68") #Macrofago
FeaturePlot(GSM3744444, features= "Cnmd") #Inmune
FeaturePlot(GSM3744444, features= "Cyp11a1") #lEY
FeaturePlot(GSM3744444, features= "Igf1") #Estromal
FeaturePlot(GSM3744444, features= "Acta2") #peritubular
GSM3744444 <-  FindNeighbors(GSM3744444, dims = 1:6)
GSM3744444 <- FindClusters(GSM3744444, resolution = c(0.3)) #11 COMUNIDADES
DimPlot(GSM3744444, reduction = "umap", label = TRUE) #RESOLUTION.0.2 == CLUSTER 10
VlnPlot(GSM3744444, features = "Ddx4")
VlnPlot(GSM3744444, features = "Dazl")
#cl_markers_GSM3744444 <- FindAllMarkers(GSM3744444,
 #                                       features = VariableFeatures(GSM3744444))
#cl_markers_GSM3744444 <- cl_markers_GSM3744444[cl_markers_GSM3744444$avg_log2FC >= 1.5 & cl_markers_GSM3744444$p_val_adj <= 0.05,]
#cat (cl_markers_GSM3744444[cl_markers_GSM3744444$cluster== "4", 7], sep = "\n")
#df_GSM3744444 <- GSM3744444@meta.data
#df = data.frame(RNA_snn_res.0.3 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
 #               cell_type = c("Sertoli cells", "Stromal cells", "Sertoli cells", "Sertoli cells", "Stromal cells", "Myoid cells", "Stromal cells", "Myoid cells",  "Germ cells", "Leydig", "Endothelial cells", "Endothelial cells"))
#df_GSM3744444$cells <- rownames(df_GSM3744444)
#df_GSM3744444 <- data.frame (merge (df_GSM3744444, df, by= "RNA_snn_res.0.3"))
#GSM3744444_UMAP <- data.frame (Embeddings(GSM3744444, reduction = "umap"))
#GSM3744444_UMAP$cells <- rownames(GSM3744444_UMAP)

#df_GSM3744444 <- merge(df_GSM3744444, GSM3744444_UMAP, by= "cells")
#write.table(df_GSM3744444, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM3744444")
#write.table(cl_markers_GSM3744444, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM3744444_markers")

GSM3744444_germinales <- subset (GSM3744444, ident= 8)
dim (GSM3744444_germinales)


############GSM4185522 - 7 días##########
GSM4185522_P7_barcodes <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM4185522_P7_Rep1/GSM4185522_P7_Rep1_barcodes.tsv.gz"
GSM4185522_P7_feature <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM4185522_P7_Rep1/GSM4185522_P7_Rep1_features.tsv.gz"
GSM4185522_P7_matrix <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM4185522_P7_Rep1/GSM4185522_P7_Rep1_matrix.mtx.gz"

#1- Creación del objeto Seurat
counts_GSM4185522 <- readMM(GSM4185522_P7_matrix)
barcodes_GSM4185522 <- read.table(GSM4185522_P7_barcodes, stringsAsFactors=F)[,1]
features_GSM4185522 <- read.csv(GSM4185522_P7_feature, stringsAsFactors=F, sep="\t", header=F)
rownames(counts_GSM4185522) <- make.unique(features_GSM4185522[,2])
colnames(counts_GSM4185522) <- barcodes_GSM4185522
GSM4185522 <- CreateSeuratObject(counts_GSM4185522, min.cells = 3, min.features = 1000, project="GSM4185522_7days")
rm (features_GSM4185522, barcodes_GSM4185522, counts_GSM4185522)
rm (GSM4185522_P7_barcodes, GSM4185522_P7_feature, GSM4185522_P7_matrix)

dim (GSM4185522)  
head (GSM4185522@meta.data)

#2- Análisis de calidad
GSM4185522[["percent.mt"]] <- PercentageFeatureSet(GSM4185522, pattern = "^mt[-\\.]")
head (GSM4185522@meta.data)
                                                   
plot1 <- FeatureScatter(GSM4185522, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(GSM4185522, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2

GSM4185522 <- subset(GSM4185522, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
plot3 <- FeatureScatter(GSM4185522, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(GSM4185522, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
rm (plot1, plot2, plot3, plot4)

dim (GSM4185522)  

#3- Eliminar multipletes (droplets)
#Creamos un archivo tipo 'Single Cell'
GSM4185522_sce <- as.SingleCellExperiment(GSM4185522)
#Corremos el scrip para la busqueda de dobletes y visualizamos el resultado
GSM4185522_sce <- scDblFinder(GSM4185522_sce, samples = "orig.ident") #Tarda unos minutos. 
head (GSM4185522_sce)
table(GSM4185522_sce$scDblFinder.class)
#Exploramos el resultado y lo incorporamos al objeto Seurat:
meta_scdblfinder_GSM4185522 <- GSM4185522_sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
head(meta_scdblfinder_GSM4185522)
rownames(meta_scdblfinder_GSM4185522) <- GSM4185522_sce@colData@rownames
GSM4185522 <- AddMetaData(object = GSM4185522, metadata = meta_scdblfinder_GSM4185522 %>% dplyr::select('scDblFinder.class'))
VlnPlot(GSM4185522, split.by = "scDblFinder.class", #ver los de sampleID
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
GSM4185522 <- subset(GSM4185522, subset = scDblFinder.class == 'singlet')
rm (meta_scdblfinder_GSM4185522, GSM4185522_sce)
GSM4185522 <- subset(GSM4185522, subset = nCount_RNA < 30000)
dim (GSM4185522)


df_GSM4185522  <- data.frame (GSM4185522@meta.data)
df_GSM4185522$QC_pass <- (df_GSM4185522$nCount_RNA < 30000 & df_GSM4185522$nFeature_RNA < 7000 & df_GSM4185522$percent.mt < 10 & df_GSM4185522$scDblFinder.class == "singlet")
write.table(df_GSM4185522, file= "~/Figuras_Paper 2026/Planillas QC/GSM4185522")

GSM4185522 <- NormalizeData(GSM4185522)
GSM4185522 <- FindVariableFeatures(GSM4185522, nfeatures = 2000)
GSM4185522 <- ScaleData(GSM4185522)
GSM4185522 <-  RunPCA(GSM4185522)
ElbowPlot(GSM4185522)

GSM4185522 <-  RunUMAP(GSM4185522, dims = 1:20) #pUEDO SEPARAR MEJOR LAS GERMINALES
UMAPPlot(GSM4185522)
FeaturePlot(GSM4185522, features= "Ddx4") #Germinal
FeaturePlot(GSM4185522, features= "Dazl") #Germinal
FeaturePlot(GSM4185522, features= "Sox9") #Sertoli
FeaturePlot(GSM4185522, features= "Amh") #Sertoli
FeaturePlot(GSM4185522, features= "Esam") #Endotelial
FeaturePlot(GSM4185522, features= "Pecam1") #Endotelial
FeaturePlot(GSM4185522, features= "Lyz2") #Macrofago
FeaturePlot(GSM4185522, features= "Cd68") #Macrofago
FeaturePlot(GSM4185522, features= "Cnmd") #Inmune
FeaturePlot(GSM4185522, features= "Cyp11a1") #lEY
FeaturePlot(GSM4185522, features= "Igf1") #Estromal
FeaturePlot(GSM4185522, features= "Acta2") #peritubular
GSM4185522 <-  FindNeighbors(GSM4185522, dims = 1:20)
GSM4185522 <- FindClusters(GSM4185522, resolution = c(1)) #11 COMUNIDADES
DimPlot(GSM4185522, reduction = "umap", label = TRUE) #RESOLUTION.0.2 == CLUSTER 10
VlnPlot(GSM4185522, features = "Ddx4")
VlnPlot(GSM4185522, features = "Dazl")
cl_markers_GSM4185522 <- FindAllMarkers(GSM4185522,
                                        features = VariableFeatures(GSM4185522))
cl_markers_GSM4185522 <- cl_markers_GSM4185522[cl_markers_GSM4185522$avg_log2FC >= 1.5 & cl_markers_GSM4185522$p_val_adj <= 0.05,]
df_GSM4185522 <- GSM4185522@meta.data
df = data.frame(RNA_snn_res.0.3 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
                cell_type = c("Sertoli cells", "Stromal cells", "Sertoli cells", "Sertoli cells", "Stromal cells", "Myoid cells", "Stromal cells", "Myoid cells",  "Germ cells", "Leydig", "Endothelial cells", "Endothelial cells"))
df_GSM4185522$cells <- rownames(df_GSM4185522)
df_GSM4185522 <- data.frame (merge (df_GSM4185522, df, by= "RNA_snn_res.0.3"))
GSM4185522_UMAP <- data.frame (Embeddings(GSM4185522, reduction = "umap"))
GSM4185522_UMAP$cells <- rownames(GSM4185522_UMAP)

df_GSM4185522 <- merge(df_GSM4185522, GSM4185522_UMAP, by= "cells")
write.table(df_GSM4185522, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM4185522")
write.table(cl_markers_GSM4185522, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM4185522_markers")

GSM4185522_germinales <- subset (GSM4185522, ident= c(6,10))
dim (GSM4185522_germinales)


#GSM4185522 <- SCTransform(GSM4185522, verbose = FALSE)

############GSM4185523 - 7 días##########
GSM4185523_P7_barcodes <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM4185523_P7_Rep2/GSM4185523_P7_Rep2_barcodes.tsv.gz"
GSM4185523_P7_feature <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM4185523_P7_Rep2/GSM4185523_P7_Rep2_features.tsv.gz"
GSM4185523_P7_matrix <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM4185523_P7_Rep2/GSM4185523_P7_Rep2_matrix.mtx.gz"

#1- Creación del objeto Seurat
counts_GSM4185523 <- readMM(GSM4185523_P7_matrix)
barcodes_GSM4185523 <- read.table(GSM4185523_P7_barcodes, stringsAsFactors=F)[,1]
features_GSM4185523 <- read.csv(GSM4185523_P7_feature, stringsAsFactors=F, sep="\t", header=F)
rownames(counts_GSM4185523) <- make.unique(features_GSM4185523[,2])
colnames(counts_GSM4185523) <- barcodes_GSM4185523
GSM4185523 <- CreateSeuratObject(counts_GSM4185523, min.cells = 3, min.features = 1000, project="GSM4185523_7days")
rm (counts_GSM4185523, barcodes_GSM4185523, features_GSM4185523)
dim (GSM4185523)  
head (GSM4185523@meta.data)
rm(GSM4185523_P7_barcodes, GSM4185523_P7_feature, GSM4185523_P7_matrix)

#2- Análisis de calidad
GSM4185523[["percent.mt"]] <- PercentageFeatureSet(GSM4185523, pattern = "^mt[-\\.]")
head (GSM4185523@meta.data)
plot5 <- FeatureScatter(GSM4185523, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot6 <- FeatureScatter(GSM4185523, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot5 + plot6

GSM4185523 <- subset(GSM4185523, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
plot7 <- FeatureScatter(GSM4185523, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot8 <- FeatureScatter(GSM4185523, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot7 + plot8

rm (plot5, plot6, plot7, plot8)
dim (GSM4185523)  

#3- Eliminar multipletes (droplets)
#Creamos un archivo tipo 'Single Cell'
GSM4185523_sce <- as.SingleCellExperiment(GSM4185523)
#Corremos el scrip para la busqueda de dobletes y visualizamos el resultado
GSM4185523_sce <- scDblFinder(GSM4185523_sce, samples = "orig.ident") #Tarda unos minutos. 
head (GSM4185523_sce)
table(GSM4185523_sce$scDblFinder.class)
#Exploramos el resultado y lo incorporamos al objeto Seurat:
meta_scdblfinder_GSM4185523 <- GSM4185523_sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
head(meta_scdblfinder_GSM4185523)
rownames(meta_scdblfinder_GSM4185523) <- GSM4185523_sce@colData@rownames
GSM4185523 <- AddMetaData(object = GSM4185523, metadata = meta_scdblfinder_GSM4185523 %>% dplyr::select('scDblFinder.class'))
VlnPlot(GSM4185523, split.by = "scDblFinder.class", #ver los de sampleID
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
GSM4185523 <- subset(GSM4185523, subset = scDblFinder.class == 'singlet')
rm (meta_scdblfinder_GSM4185523, GSM4185523_sce)
GSM4185523 <- subset(GSM4185523, subset = nCount_RNA < 30000)
dim (GSM4185523)





df_GSM4185523  <- data.frame (GSM4185523@meta.data)
df_GSM4185523$QC_pass <- (df_GSM4185523$nCount_RNA < 30000 & df_GSM4185523$nFeature_RNA < 7000 & df_GSM4185523$percent.mt < 10 & df_GSM4185523$scDblFinder.class == "singlet")
write.table(df_GSM4185523, file= "~/Figuras_Paper 2026/Planillas QC/GSM4185523")

GSM4185523 <- NormalizeData(GSM4185523)
GSM4185523 <- FindVariableFeatures(GSM4185523, nfeatures = 2000)
GSM4185523 <- ScaleData(GSM4185523)
GSM4185523 <-  RunPCA(GSM4185523)
ElbowPlot(GSM4185523)

GSM4185523 <-  RunUMAP(GSM4185523, dims = 1:20) #pUEDO SEPARAR MEJOR LAS GERMINALES
UMAPPlot(GSM4185523)
FeaturePlot(GSM4185523, features= "Ddx4") #Germinal
FeaturePlot(GSM4185523, features= "Dazl") #Germinal
FeaturePlot(GSM4185523, features= "Sox9") #Sertoli
FeaturePlot(GSM4185523, features= "Amh") #Sertoli
FeaturePlot(GSM4185523, features= "Esam") #Endotelial
FeaturePlot(GSM4185523, features= "Pecam1") #Endotelial
FeaturePlot(GSM4185523, features= "Lyz2") #Macrofago
FeaturePlot(GSM4185523, features= "Cd68") #Macrofago
FeaturePlot(GSM4185523, features= "Cnmd") #Inmune
FeaturePlot(GSM4185523, features= "Cyp11a1") #lEY
FeaturePlot(GSM4185523, features= "Igf1") #Estromal
FeaturePlot(GSM4185523, features= "Acta2") #peritubular
GSM4185523 <-  FindNeighbors(GSM4185523, dims = 1:20)
GSM4185523 <- FindClusters(GSM4185523, resolution = c(0.3)) #11 COMUNIDADES
DimPlot(GSM4185523, reduction = "umap", label = TRUE) #RESOLUTION.0.2 == CLUSTER 10
VlnPlot(GSM4185523, features = "Ddx4")
VlnPlot(GSM4185523, features = "Dazl")
cl_markers_GSM4185523 <- FindAllMarkers(GSM4185523,
                                        features = VariableFeatures(GSM4185523))
cl_markers_GSM4185523 <- cl_markers_GSM4185523[cl_markers_GSM4185523$avg_log2FC >= 1.5 & cl_markers_GSM4185523$p_val_adj <= 0.05,]
df_GSM4185523 <- GSM4185523@meta.data
df = data.frame(RNA_snn_res.0.3 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
                cell_type = c("Sertoli cells", "Stromal cells", "Sertoli cells", "Sertoli cells", "Stromal cells", "Myoid cells", "Stromal cells", "Myoid cells",  "Germ cells", "Leydig", "Endothelial cells", "Endothelial cells"))
df_GSM4185523$cells <- rownames(df_GSM4185523)
df_GSM4185523 <- data.frame (merge (df_GSM4185523, df, by= "RNA_snn_res.0.3"))
GSM4185523_UMAP <- data.frame (Embeddings(GSM4185523, reduction = "umap"))
GSM4185523_UMAP$cells <- rownames(GSM4185523_UMAP)

df_GSM4185523 <- merge(df_GSM4185523, GSM4185523_UMAP, by= "cells")
write.table(df_GSM4185523, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM4185523")
write.table(cl_markers_GSM4185523, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM4185523_markers")

GSM4185523_germinales <- subset (GSM4185523, ident= c(5,7,8))
dim (GSM4185523_germinales)
#GSM4185523 <- SCTransform(GSM4185523, verbose = FALSE)




############ GSM5708025 - 14 días##########
GSM5708025_P14_barcodes <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM5708025_P14/GSM5708025_vivo_14_barcodes.tsv.gz"
GSM5708025_P14_feature <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM5708025_P14/GSM5708025_vivo_14_features.tsv.gz"
GSM5708025_P14_matrix <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM5708025_P14/GSM5708025_vivo_14_matrix.mtx.gz"  

#1- Creación del objeto Seurat
counts_GSM5708025 <- readMM(GSM5708025_P14_matrix)
barcodes_GSM5708025 <- read.table(GSM5708025_P14_barcodes, stringsAsFactors=F)[,1]
features_GSM5708025 <- read.csv(GSM5708025_P14_feature, stringsAsFactors=F, sep="\t", header=F)
rownames(counts_GSM5708025) <- make.unique(features_GSM5708025[,2])
colnames(counts_GSM5708025) <- barcodes_GSM5708025
GSM5708025 <- CreateSeuratObject(counts_GSM5708025, min.cells = 3, min.features = 1000, project="GSM5708025_14days")
rm (features_GSM5708025, barcodes_GSM5708025, counts_GSM5708025)
rm(GSM5708025_P14_barcodes, GSM5708025_P14_feature, GSM5708025_P14_matrix)

dim (GSM5708025)  
head (GSM5708025@meta.data)

#2- Análisis de calidad
GSM5708025[["percent.mt"]] <- PercentageFeatureSet(GSM5708025, pattern = "^mt[-\\.]")
head (GSM5708025@meta.data)

plot1 <- FeatureScatter(GSM5708025, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(GSM5708025, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2

GSM5708025 <- subset(GSM5708025, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
plot3 <- FeatureScatter(GSM5708025, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(GSM5708025, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
rm (plot1, plot2, plot3, plot4)

dim (GSM5708025)  

#3- Eliminar multipletes (droplets)
#Creamos un archivo tipo 'Single Cell'
GSM5708025_sce <- as.SingleCellExperiment(GSM5708025)
#Corremos el scrip para la busqueda de dobletes y visualizamos el resultado
GSM5708025_sce <- scDblFinder(GSM5708025_sce, samples = "orig.ident") #Tarda unos minutos. 
head (GSM5708025_sce)
table(GSM5708025_sce$scDblFinder.class)
#Exploramos el resultado y lo incorporamos al objeto Seurat:
meta_scdblfinder_GSM5708025 <- GSM5708025_sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
head(meta_scdblfinder_GSM5708025)
rownames(meta_scdblfinder_GSM5708025) <- GSM5708025_sce@colData@rownames
GSM5708025 <- AddMetaData(object = GSM5708025, metadata = meta_scdblfinder_GSM5708025 %>% dplyr::select('scDblFinder.class'))
VlnPlot(GSM5708025, split.by = "scDblFinder.class", #ver los de sampleID
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
GSM5708025 <- subset(GSM5708025, subset = scDblFinder.class == 'singlet')
rm (meta_scdblfinder_GSM5708025, GSM5708025_sce)
GSM5708025 <- subset(GSM5708025, subset = nCount_RNA < 30000)
dim (GSM5708025)

df_GSM5708025  <- data.frame (GSM5708025@meta.data)
df_GSM5708025$QC_pass <- (df_GSM5708025$nCount_RNA < 30000 & df_GSM5708025$nFeature_RNA < 7000 & df_GSM5708025$percent.mt < 10 & df_GSM5708025$scDblFinder.class == "singlet")
write.table(df_GSM5708025, file= "~/Figuras_Paper 2026/Planillas QC/GSM5708025")

GSM5708025 <- NormalizeData(GSM5708025)
GSM5708025 <- FindVariableFeatures(GSM5708025, nfeatures = 2000)
GSM5708025 <- ScaleData(GSM5708025)
GSM5708025 <-  RunPCA(GSM5708025)
ElbowPlot(GSM5708025)

GSM5708025 <-  RunUMAP(GSM5708025, dims = 1:12) #pUEDO SEPARAR MEJOR LAS GERMINALES
UMAPPlot(GSM5708025)
FeaturePlot(GSM5708025, features= "Ddx4") #Germinal
FeaturePlot(GSM5708025, features= "Dazl") #Germinal
FeaturePlot(GSM5708025, features= "Sox9") #Sertoli
FeaturePlot(GSM5708025, features= "Amh") #Sertoli
FeaturePlot(GSM5708025, features= "Esam") #Endotelial
FeaturePlot(GSM5708025, features= "Pecam1") #Endotelial
FeaturePlot(GSM5708025, features= "Lyz2") #Macrofago
FeaturePlot(GSM5708025, features= "Cd68") #Macrofago
FeaturePlot(GSM5708025, features= "Cnmd") #Inmune
FeaturePlot(GSM5708025, features= "Cyp11a1") #lEY
FeaturePlot(GSM5708025, features= "Igf1") #Estromal
FeaturePlot(GSM5708025, features= "Acta2") #peritubular
GSM5708025 <-  FindNeighbors(GSM5708025, dims = 1:12)
GSM5708025 <- FindClusters(GSM5708025, resolution = c(0.3)) #11 COMUNIDADES
DimPlot(GSM5708025, reduction = "umap", label = TRUE) #RESOLUTION.0.2 == CLUSTER 10
VlnPlot(GSM5708025, features = "Ddx4")
VlnPlot(GSM5708025, features = "Dazl")
cl_markers_GSM5708025 <- FindAllMarkers(GSM5708025,
                                        features = VariableFeatures(GSM5708025))
cl_markers_GSM5708025 <- cl_markers_GSM5708025[cl_markers_GSM5708025$avg_log2FC >= 1.5 & cl_markers_GSM5708025$p_val_adj <= 0.05,]
df_GSM5708025 <- GSM5708025@meta.data
df = data.frame(RNA_snn_res.0.3 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
                cell_type = c("Sertoli cells", "Stromal cells", "Sertoli cells", "Sertoli cells", "Stromal cells", "Myoid cells", "Stromal cells", "Myoid cells",  "Germ cells", "Leydig", "Endothelial cells", "Endothelial cells"))
df_GSM5708025$cells <- rownames(df_GSM5708025)
df_GSM5708025 <- data.frame (merge (df_GSM5708025, df, by= "RNA_snn_res.0.3"))
GSM5708025_UMAP <- data.frame (Embeddings(GSM5708025, reduction = "umap"))
GSM5708025_UMAP$cells <- rownames(GSM5708025_UMAP)

df_GSM5708025 <- merge(df_GSM5708025, GSM5708025_UMAP, by= "cells")
write.table(df_GSM5708025, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM5708025")
write.table(cl_markers_GSM5708025, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM5708025_markers")

GSM5708025_germinales <- subset (GSM5708025, ident= c(0,1,2,5,6,8))
dim (GSM5708025_germinales)

#GSM5708025 <- SCTransform(GSM5708025, verbose = FALSE)

############GSM6468366 - 21 días##########
GSM6468366_P21 <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM6468366_P21/GSM6468366_WT-Day21.raw_feature_bc_matrix.h5"

#1- Creación del objeto Seurat
GSM6468366 <- Read10X_h5(GSM6468366_P21)
GSM6468366 <- CreateSeuratObject(counts = GSM6468366, min.cells = 3, min.features = 1000, project = "GSM6468366_21days")
dim (GSM6468366) #32286 2135905
head (GSM6468366@meta.data)
rm (GSM6468366_P21)

#2- Análisis de calidad
GSM6468366[["percent.mt"]] <- PercentageFeatureSet(GSM6468366, pattern = "^mt[-\\.]")
head (GSM6468366@meta.data)
plot13 <- FeatureScatter(GSM6468366, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot14 <- FeatureScatter(GSM6468366, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot13 + plot14

GSM6468366 <- subset(GSM6468366, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10) #ACA DEJE MAS ADN MITOCODRIAL- VER SI NO AFECTA
dim (GSM6468366) #32286  2832
plot15 <- FeatureScatter(GSM6468366, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot16 <- FeatureScatter(GSM6468366, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot15 + plot16
rm (plot13, plot14, plot15, plot16)

#3- Eliminar multipletes (droplets)
#Creamos un archivo tipo 'Single Cell'
GSM6468366_sce <- as.SingleCellExperiment(GSM6468366)
#Corremos el scrip para la busqueda de dobletes y visualizamos el resultado
GSM6468366_sce <- scDblFinder(GSM6468366_sce, samples = "orig.ident") #Tarda unos minutos. 
head (GSM6468366_sce)
table(GSM6468366_sce$scDblFinder.class)
#Exploramos el resultado y lo incorporamos al objeto Seurat:
meta_scdblfinder_GSM6468366 <- GSM6468366_sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
head(meta_scdblfinder_GSM6468366)
rownames(meta_scdblfinder_GSM6468366) <- GSM6468366_sce@colData@rownames
GSM6468366 <- AddMetaData(object = GSM6468366, metadata = meta_scdblfinder_GSM6468366 %>% dplyr::select('scDblFinder.class'))
VlnPlot(GSM6468366, split.by = "scDblFinder.class", #ver los de sampleID
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
GSM6468366 <- subset(GSM6468366, subset = scDblFinder.class == 'singlet')
rm (meta_scdblfinder_GSM6468366, GSM6468366_sce)


GSM6468366 <- subset(GSM6468366, subset = nCount_RNA < 30000)
dim (GSM6468366)
df_GSM6468366  <- data.frame (GSM6468366@meta.data)
df_GSM6468366$QC_pass <- (df_GSM6468366$nCount_RNA < 30000 & df_GSM6468366$nFeature_RNA < 7000 & df_GSM6468366$percent.mt < 10 & df_GSM6468366$scDblFinder.class == "singlet")
write.table(df_GSM6468366, file= "~/Figuras_Paper 2026/Planillas QC/GSM6468366")

GSM6468366 <- NormalizeData(GSM6468366)
GSM6468366 <- FindVariableFeatures(GSM6468366, nfeatures = 2000)
GSM6468366 <- ScaleData(GSM6468366)
GSM6468366 <-  RunPCA(GSM6468366)
ElbowPlot(GSM6468366)

GSM6468366 <-  RunUMAP(GSM6468366, dims = 1:5) #pUEDO SEPARAR MEJOR LAS GERMINALES
UMAPPlot(GSM6468366)
FeaturePlot(GSM6468366, features= "Ddx4") #Germinal
FeaturePlot(GSM6468366, features= "Dazl") #Germinal
FeaturePlot(GSM6468366, features= "Sycp3") #Germinal
FeaturePlot(GSM6468366, features= "Sox9") #Sertoli
FeaturePlot(GSM6468366, features= "Amh") #Sertoli
FeaturePlot(GSM6468366, features= "Esam") #Endotelial
FeaturePlot(GSM6468366, features= "Pecam1") #Endotelial
FeaturePlot(GSM6468366, features= "Lyz2") #Macrofago
FeaturePlot(GSM6468366, features= "Cd68") #Macrofago
FeaturePlot(GSM6468366, features= "Cnmd") #Inmune
FeaturePlot(GSM6468366, features= "Cyp11a1") #lEY
FeaturePlot(GSM6468366, features= "Igf1") #Estromal
FeaturePlot(GSM6468366, features= "Acta2") #peritubular
GSM6468366 <-  FindNeighbors(GSM6468366, dims = 1:5)
GSM6468366 <- FindClusters(GSM6468366, resolution = c(0.3)) #11 COMUNIDADES
DimPlot(GSM6468366, reduction = "umap", label = TRUE) #RESOLUTION.0.2 == CLUSTER 10
VlnPlot(GSM6468366, features = "Ddx4")
VlnPlot(GSM6468366, features = "Sox9")
cl_markers_GSM6468366 <- FindAllMarkers(GSM6468366,
                                        features = VariableFeatures(GSM6468366))
cl_markers_GSM6468366 <- cl_markers_GSM6468366[cl_markers_GSM6468366$avg_log2FC >= 1.5 & cl_markers_GSM6468366$p_val_adj <= 0.05,]
df_GSM6468366 <- GSM6468366@meta.data
df = data.frame(RNA_snn_res.0.3 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
                cell_type = c("Sertoli cells", "Stromal cells", "Sertoli cells", "Sertoli cells", "Stromal cells", "Myoid cells", "Stromal cells", "Myoid cells",  "Germ cells", "Leydig", "Endothelial cells", "Endothelial cells"))
df_GSM6468366$cells <- rownames(df_GSM6468366)
df_GSM6468366 <- data.frame (merge (df_GSM6468366, df, by= "RNA_snn_res.0.3"))
GSM6468366_UMAP <- data.frame (Embeddings(GSM6468366, reduction = "umap"))
GSM6468366_UMAP$cells <- rownames(GSM6468366_UMAP)

df_GSM6468366 <- merge(df_GSM6468366, GSM6468366_UMAP, by= "cells")
write.table(df_GSM6468366, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM6468366")
write.table(cl_markers_GSM6468366, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM6468366_markers")

GSM6468366_germinales <- subset (GSM6468366, ident= c(0,1,2,3,6))
dim (GSM6468366_germinales)


############GSM7135562 - 24días ##########
GSM7135562_barcodes <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM7135562_P24/GSM7135562_barcodes.tsv.gz"
GSM7135562_feature <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM7135562_P24/GSM7135562_features.tsv.gz"
GSM7135562_matrix <- "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/Datos/GSM7135562_P24/GSM7135562_matrix.mtx.gz"  

#1- Creación del objeto Seurat
counts_GSM7135562 <- readMM(GSM7135562_matrix)
barcodes_GSM7135562 <- read.table(GSM7135562_barcodes, stringsAsFactors=F)[,1]
features_GSM7135562 <- read.csv(GSM7135562_feature, stringsAsFactors=F, sep="\t", header=F)
rownames(counts_GSM7135562) <- make.unique(features_GSM7135562[,2])
colnames(counts_GSM7135562) <- barcodes_GSM7135562
GSM7135562 <- CreateSeuratObject(counts_GSM7135562, min.cells = 3, min.features = 1000, project="GSM7135562_24days")
rm (GSM7135562_barcodes, GSM7135562_feature, GSM7135562_matrix, counts_GSM7135562, barcodes_GSM7135562, features_GSM7135562)

dim (GSM7135562)  
head (GSM7135562@meta.data)

#2- Análisis de calidad
GSM7135562[["percent.mt"]] <- PercentageFeatureSet(GSM7135562, pattern = "^mt[-\\.]")
head (GSM7135562@meta.data)

plot1 <- FeatureScatter(GSM7135562, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(GSM7135562, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2

GSM7135562 <- subset(GSM7135562, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10 & nCount_RNA <30000)
plot3 <- FeatureScatter(GSM7135562, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(GSM7135562, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
rm (plot1, plot2, plot3, plot4)

dim (GSM7135562)  #23327  2450

#3- Eliminar multipletes (droplets)
#Creamos un archivo tipo 'Single Cell'
GSM7135562_sce <- as.SingleCellExperiment(GSM7135562)
#Corremos el scrip para la busqueda de dobletes y visualizamos el resultado
GSM7135562_sce <- scDblFinder(GSM7135562_sce, samples = "orig.ident") #Tarda unos minutos. 
head (GSM7135562_sce)
table(GSM7135562_sce$scDblFinder.class) #71 droplets
#Exploramos el resultado y lo incorporamos al objeto Seurat:
meta_scdblfinder_GSM7135562 <- GSM7135562_sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
head(meta_scdblfinder_GSM7135562)
rownames(meta_scdblfinder_GSM7135562) <- GSM7135562_sce@colData@rownames
GSM7135562 <- AddMetaData(object = GSM7135562, metadata = meta_scdblfinder_GSM7135562 %>% dplyr::select('scDblFinder.class'))
VlnPlot(GSM7135562, split.by = "scDblFinder.class", #ver los de sampleID
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
GSM7135562 <- subset(GSM7135562, subset = scDblFinder.class == 'singlet')
rm (meta_scdblfinder_GSM7135562, GSM7135562_sce)
GSM7135562 <- subset(GSM7135562, subset = nCount_RNA < 30000)
dim (GSM7135562) 

df_GSM7135562  <- data.frame (GSM7135562@meta.data)
df_GSM7135562$QC_pass <- (df_GSM7135562$nCount_RNA < 30000 & df_GSM7135562$nFeature_RNA < 7000 & df_GSM7135562$percent.mt < 10 & df_GSM7135562$scDblFinder.class == "singlet")
write.table(df_GSM7135562, file= "~/Figuras_Paper 2026/Planillas QC/GSM7135562")

GSM7135562 <- NormalizeData(GSM7135562)
GSM7135562 <- FindVariableFeatures(GSM7135562, nfeatures = 2000)
GSM7135562 <- ScaleData(GSM7135562)
GSM7135562 <-  RunPCA(GSM7135562)
ElbowPlot(GSM7135562)

GSM7135562 <-  RunUMAP(GSM7135562, dims = 1:20) #pUEDO SEPARAR MEJOR LAS GERMINALES
UMAPPlot(GSM7135562)
FeaturePlot(GSM7135562, features= "Ddx4") #Germinal
FeaturePlot(GSM7135562, features= "Dazl") #Germinal
FeaturePlot(GSM7135562, features= "Sox9") #Sertoli
FeaturePlot(GSM7135562, features= "Amh") #Sertoli
FeaturePlot(GSM7135562, features= "Esam") #Endotelial
FeaturePlot(GSM7135562, features= "Pecam1") #Endotelial
FeaturePlot(GSM7135562, features= "Lyz2") #Macrofago
FeaturePlot(GSM7135562, features= "Cd68") #Macrofago
FeaturePlot(GSM7135562, features= "Cnmd") #Inmune
FeaturePlot(GSM7135562, features= "Cyp11a1") #lEY
FeaturePlot(GSM7135562, features= "Igf1") #Estromal
FeaturePlot(GSM7135562, features= "Acta2") #peritubular
GSM7135562 <-  FindNeighbors(GSM7135562, dims = 1:3)
GSM7135562 <- FindClusters(GSM7135562, resolution = c(0.1)) #11 COMUNIDADES
DimPlot(GSM7135562, reduction = "umap", label = TRUE) #RESOLUTION.0.2 == CLUSTER 10
VlnPlot(GSM7135562, features = "Ddx4")
VlnPlot(GSM7135562, features = "Dazl")
cl_markers_GSM7135562 <- FindAllMarkers(GSM7135562,
                                        features = VariableFeatures(GSM7135562))
cl_markers_GSM7135562 <- cl_markers_GSM7135562[cl_markers_GSM7135562$avg_log2FC >= 1.5 & cl_markers_GSM7135562$p_val_adj <= 0.05,]
df_GSM7135562 <- GSM7135562@meta.data
df = data.frame(RNA_snn_res.0.3 = c("0", "1", "2", "3", "4", "5", "6", "8", "9", "10", "11"), 
                cell_type = c("Sertoli cells", "Stromal cells", "Sertoli cells", "Sertoli cells", "Stromal cells", "Myoid cells", "Stromal cells", "Myoid cells",  "Germ cells", "Leydig", "Endothelial cells", "Endothelial cells"))
df_GSM7135562$cells <- rownames(df_GSM7135562)
df_GSM7135562 <- data.frame (merge (df_GSM7135562, df, by= "RNA_snn_res.0.3"))
GSM7135562_UMAP <- data.frame (Embeddings(GSM7135562, reduction = "umap"))
GSM7135562_UMAP$cells <- rownames(GSM7135562_UMAP)

df_GSM7135562 <- merge(df_GSM7135562, GSM7135562_UMAP, by= "cells")
write.table(df_GSM7135562, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM7135562")
write.table(cl_markers_GSM7135562, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSM7135562_markers")

GSM7135562_germinales <- GSM7135562
dim (GSM7135562_germinales)



###############SRR6129051-Adulto 2###############
SRR6129051_A_barcodes <- "~/SRR6129051_filtered_feature_bc_matrix/barcodes.tsv.gz"
SRR6129051_A_feature <- "~/SRR6129051_filtered_feature_bc_matrix/features.tsv.gz"
SRR6129051_A_matrix <- "~/SRR6129051_filtered_feature_bc_matrix//matrix.mtx.gz"  

#1- Creación del objeto Seurat
counts_SRR6129051 <- readMM(SRR6129051_A_matrix)
barcodes_SRR6129051 <- read.table(SRR6129051_A_barcodes, stringsAsFactors=F)[,1]
features_SRR6129051 <- read.csv(SRR6129051_A_feature, stringsAsFactors=F, sep="\t", header=F)
rownames(counts_SRR6129051) <- make.unique(features_SRR6129051[,2])
colnames(counts_SRR6129051) <- barcodes_SRR6129051
SRR6129051 <- CreateSeuratObject(counts_SRR6129051, min.cells = 3, min.features = 1000, project="SRR6129051_8weeks")
rm (SRR6129051_A_barcodes, SRR6129051_A_feature, SRR6129051_A_matrix, counts_SRR6129051, barcodes_SRR6129051, features_SRR6129051)
dim (SRR6129051)  

#2- Análisis de calidad
SRR6129051[["percent.mt"]] <- PercentageFeatureSet(SRR6129051, pattern = "^mt[-\\.]")
head (SRR6129051@meta.data)
summary (SRR6129051@meta.data)

plot1 <- FeatureScatter(SRR6129051, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(SRR6129051, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2

SRR6129051 <- subset(SRR6129051, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
dim (SRR6129051)

#Creamos un archivo tipo 'Single Cell'
SRR6129051_sce <- as.SingleCellExperiment(SRR6129051)
#Corremos el scrip para la busqueda de dobletes y visualizamos el resultado
SRR6129051_sce <- scDblFinder(SRR6129051_sce, samples = "orig.ident") #Tarda unos minutos. 
table(SRR6129051_sce$scDblFinder.class) #71 droplets
#Exploramos el resultado y lo incorporamos al objeto Seurat:
meta_scdblfinder_SRR6129051 <- SRR6129051_sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
head(meta_scdblfinder_SRR6129051)
rownames(meta_scdblfinder_SRR6129051) <- SRR6129051_sce@colData@rownames
SRR6129051 <- AddMetaData(object = SRR6129051, metadata = meta_scdblfinder_SRR6129051 %>% dplyr::select('scDblFinder.class'))
VlnPlot(SRR6129051, split.by = "scDblFinder.class", #ver los de sampleID
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
SRR6129051 <- subset(SRR6129051 , subset = scDblFinder.class == 'singlet')
SRR6129051 <- subset(SRR6129051, subset = nCount_RNA < 30000)
dim (SRR6129051) #[1] 22820  4115

df_SRR6129051  <- data.frame (SRR6129051@meta.data)
df_SRR6129051$QC_pass <- (df_SRR6129051$nCount_RNA < 30000 & df_SRR6129051$nFeature_RNA < 7000 & df_SRR6129051$percent.mt < 10 & df_SRR6129051$scDblFinder.class == "singlet")
write.table(df_SRR6129051, file= "~/Figuras_Paper 2026/Planillas QC/SRR6129051")

SRR6129051 <- NormalizeData(SRR6129051)
SRR6129051 <- FindVariableFeatures(SRR6129051, nfeatures = 2000)
SRR6129051 <- ScaleData(SRR6129051)
SRR6129051 <-  RunPCA(SRR6129051)
ElbowPlot(SRR6129051)

SRR6129051 <-  RunUMAP(SRR6129051, dims = 1:5) #pUEDO SEPARAR MEJOR LAS GERMINALES
UMAPPlot(SRR6129051)
FeaturePlot(SRR6129051, features= "Ddx4") #Germinal
FeaturePlot(SRR6129051, features= "Dazl") #Germinal
FeaturePlot(SRR6129051, features= "Sox9") #Sertoli
FeaturePlot(SRR6129051, features= "Amh") #Sertoli
FeaturePlot(SRR6129051, features= "Esam") #Endotelial
FeaturePlot(SRR6129051, features= "Pecam1") #Endotelial
FeaturePlot(SRR6129051, features= "Lyz2") #Macrofago
FeaturePlot(SRR6129051, features= "Cd68") #Macrofago
FeaturePlot(SRR6129051, features= "Cnmd") #Inmune
FeaturePlot(SRR6129051, features= "Cyp11a1") #lEY
FeaturePlot(SRR6129051, features= "Igf1") #Estromal
FeaturePlot(SRR6129051, features= "nFeature_RNA") #peritubular
SRR6129051 <-  FindNeighbors(SRR6129051, dims = 1:3)
SRR6129051 <- FindClusters(SRR6129051, resolution = c(0.3)) #11 COMUNIDADES
DimPlot(SRR6129051, reduction = "umap", label = TRUE) #RESOLUTION.0.2 == CLUSTER 10
FeaturePlot(SRR6129051, features = "nFeature_RNA")
VlnPlot(SRR6129051, features = "Dazl")
cl_markers_SRR6129051 <- FindAllMarkers(SRR6129051,
                                        features = VariableFeatures(SRR6129051))
cl_markers_SRR6129051 <- cl_markers_SRR6129051[cl_markers_SRR6129051$avg_log2FC >= 1.5 & cl_markers_SRR6129051$p_val_adj <= 0.05,]
df_SRR6129051 <- SRR6129051@meta.data
df = data.frame(RNA_snn_res.0.3 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
                cell_type = c("Sertoli cells", "Stromal cells", "Sertoli cells", "Sertoli cells", "Stromal cells", "Myoid cells", "Stromal cells", "Myoid cells",  "Germ cells", "Leydig", "Endothelial cells", "Endothelial cells"))
df_SRR6129051$cells <- rownames(df_SRR6129051)
df_SRR6129051 <- data.frame (merge (df_SRR6129051, df, by= "RNA_snn_res.0.3"))
SRR6129051_UMAP <- data.frame (Embeddings(SRR6129051, reduction = "umap"))
SRR6129051_UMAP$cells <- rownames(SRR6129051_UMAP)

df_SRR6129051 <- merge(df_SRR6129051, SRR6129051_UMAP, by= "cells")
write.table(df_SRR6129051, file= "~/Figuras_Paper 2026/UMAP y Clustering/SRR6129051")
write.table(cl_markers_SRR6129051, file= "~/Figuras_Paper 2026/UMAP y Clustering/SRR6129051_markers")

SRR6129051_germinales <- SRR6129051
dim (SRR6129051_germinales)





###############SRR6129050-Adulto 2###############
SRR6129050_A_barcodes <- "~/SRR6129050_filtered_feature_bc_matrix/barcodes.tsv.gz"
SRR6129050_A_feature <- "~/SRR6129050_filtered_feature_bc_matrix/features.tsv.gz"
SRR6129050_A_matrix <- "~/SRR6129050_filtered_feature_bc_matrix//matrix.mtx.gz"  

#1- Creación del objeto Seurat
counts_SRR6129050 <- readMM(SRR6129050_A_matrix)
barcodes_SRR6129050 <- read.table(SRR6129050_A_barcodes, stringsAsFactors=F)[,1]
features_SRR6129050 <- read.csv(SRR6129050_A_feature, stringsAsFactors=F, sep="\t", header=F)
rownames(counts_SRR6129050) <- make.unique(features_SRR6129050[,2])
colnames(counts_SRR6129050) <- barcodes_SRR6129050
SRR6129050 <- CreateSeuratObject(counts_SRR6129050, min.cells = 3, min.features = 1000, project="SRR6129050_8weeks")
rm (SRR6129050_A_barcodes, SRR6129050_A_feature, SRR6129050_A_matrix, counts_SRR6129050, barcodes_SRR6129050, features_SRR6129050)

dim (SRR6129050)  

#2- Análisis de calidad
SRR6129050[["percent.mt"]] <- PercentageFeatureSet(SRR6129050, pattern = "^mt[-\\.]")
head (SRR6129050@meta.data)
summary (SRR6129050@meta.data)

plot1 <- FeatureScatter(SRR6129050, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(SRR6129050, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2

SRR6129050<- subset(SRR6129050, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
dim (SRR6129050)

#Creamos un archivo tipo 'Single Cell'
SRR6129050_sce <- as.SingleCellExperiment(SRR6129050)
#Corremos el scrip para la busqueda de dobletes y visualizamos el resultado
SRR6129050_sce <- scDblFinder(SRR6129050_sce, samples = "orig.ident") #Tarda unos minutos. 
table(SRR6129050_sce$scDblFinder.class) #71 droplets
#Exploramos el resultado y lo incorporamos al objeto Seurat:
meta_scdblfinder_SRR6129050 <- SRR6129050_sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
head(meta_scdblfinder_SRR6129050)
rownames(meta_scdblfinder_SRR6129050) <- SRR6129050_sce@colData@rownames
SRR6129050 <- AddMetaData(object = SRR6129050, metadata = meta_scdblfinder_SRR6129050 %>% dplyr::select('scDblFinder.class'))
VlnPlot(SRR6129050, split.by = "scDblFinder.class", #ver los de sampleID
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
SRR6129050<- subset(SRR6129050, subset = scDblFinder.class == 'singlet')
dim (SRR6129050) #[1] 22820  4115
SRR6129050 <- subset(SRR6129050, subset = nCount_RNA < 30000)

SRR6129050 <- NormalizeData(SRR6129050)
SRR6129050 <- FindVariableFeatures(SRR6129050, nfeatures = 2000)
SRR6129050 <- ScaleData(SRR6129050)
SRR6129050 <-  RunPCA(SRR6129050)
ElbowPlot(SRR6129050)

SRR6129050 <-  RunUMAP(SRR6129050, dims = 1:3) #pUEDO SEPARAR MEJOR LAS GERMINALES
UMAPPlot(SRR6129050)
FeaturePlot(SRR6129050, features= "Ddx4") #Germinal
FeaturePlot(SRR6129050, features= "nFeature_RNA") #Germinal

SRR6129050 <-  FindNeighbors(SRR6129050, dims = 1:3)
SRR6129050 <- FindClusters(SRR6129050, resolution = c(0.3)) #11 COMUNIDADES
DimPlot(SRR6129050, reduction = "umap", label = TRUE) #RESOLUTION.0.2 == CLUSTER 10

VlnPlot(SRR6129050, features = "Ddx4")

SRR6129050_germinales <- SRR6129050
dim (SRR6129050_germinales)


####################GSE221226##############################
GSE221226_barcodes <- "C:/Users/Usuario/Downloads/GSM6856439_WT_barcodes.tsv.gz"
GSE221226_feature <- "C:/Users/Usuario/Downloads/GSM6856439_WT_genes.tsv.gz"
GSE221226_matrix <- "C:/Users/Usuario/Downloads/GSM6856439_WT_matrix.mtx.gz"  

#1- Creación del objeto Seurat
counts_GSE221226 <- readMM(GSE221226_matrix)
barcodes_GSE221226 <- read.table(GSE221226_barcodes, stringsAsFactors=F)[,1]
features_GSE221226 <- read.csv(GSE221226_feature, stringsAsFactors=F, sep="\t", header=F)
rownames(counts_GSE221226) <- make.unique(features_GSE221226[,2])
colnames(counts_GSE221226) <- barcodes_GSE221226
GSE221226 <- CreateSeuratObject(counts_GSE221226, min.cells = 3, min.features = 1000, project="GSE221226_8weeks")
rm (GSE221226_barcodes, GSE221226_feature, GSE221226_matrix, counts_GSE221226, barcodes_GSE221226, features_GSE221226)

dim (GSE221226)  
head (GSE221226@meta.data)

#2- Análisis de calidad
GSE221226[["percent.mt"]] <- PercentageFeatureSet(GSE221226, pattern = "^mt[-\\.]")
head (GSE221226@meta.data)

plot1 <- FeatureScatter(GSE221226, feature1 = "nCount_RNA", feature2 = "percent.mt") 
plot2 <- FeatureScatter(GSE221226, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
plot1 + plot2

GSE221226 <- subset(GSE221226, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
plot3 <- FeatureScatter(GSE221226, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot4 <- FeatureScatter(GSE221226, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 + plot4
rm (plot1, plot2, plot3, plot4)

dim (GSE221226)  #23327  2450

#3- Eliminar multipletes (droplets)
#Creamos un archivo tipo 'Single Cell'
GSE221226_sce <- as.SingleCellExperiment(GSE221226)
#Corremos el scrip para la busqueda de dobletes y visualizamos el resultado
GSE221226_sce <- scDblFinder(GSE221226_sce, samples = "orig.ident") #Tarda unos minutos. 
head (GSE221226_sce)
table(GSE221226_sce$scDblFinder.class) #71 droplets
#Exploramos el resultado y lo incorporamos al objeto Seurat:
meta_scdblfinder_GSE221226 <- GSE221226_sce@colData@listData %>% as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
head(meta_scdblfinder_GSE221226)
rownames(meta_scdblfinder_GSE221226) <- GSE221226_sce@colData@rownames
GSE221226 <- AddMetaData(object = GSE221226, metadata = meta_scdblfinder_GSE221226 %>% dplyr::select('scDblFinder.class'))
VlnPlot(GSE221226, split.by = "scDblFinder.class", #ver los de sampleID
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0) + theme(legend.position = 'right')
GSE221226 <- subset(GSE221226, subset = scDblFinder.class == 'singlet')
rm (meta_scdblfinder_GSE221226, GSE221226_sce)
GSE221226 <- subset(GSE221226, subset = nCount_RNA < 30000)
dim (GSE221226) #[1] 23327  2379

df_GSE221226  <- data.frame (GSE221226@meta.data)
df_GSE221226$QC_pass <- (df_GSE221226$nCount_RNA < 30000 & df_GSE221226$nFeature_RNA < 7000 & df_GSE221226$percent.mt < 10 & df_GSE221226$scDblFinder.class == "singlet")
write.table(df_GSE221226, file= "~/Figuras_Paper 2026/Planillas QC/GSE221226")

GSE221226 <- NormalizeData(GSE221226)
GSE221226 <- FindVariableFeatures(GSE221226, nfeatures = 2000)
GSE221226 <- ScaleData(GSE221226)
GSE221226 <-  RunPCA(GSE221226)
ElbowPlot(GSE221226)

GSE221226 <-  RunUMAP(GSE221226, dims = 1:5) #pUEDO SEPARAR MEJOR LAS GERMINALES
UMAPPlot(GSE221226)
FeaturePlot(GSE221226, features= "Ddx4") #Germinal
FeaturePlot(GSE221226, features= "Dazl") #Germinal
FeaturePlot(GSE221226, features= "Sox9") #Sertoli
FeaturePlot(GSE221226, features= "Amh") #Sertoli
FeaturePlot(GSE221226, features= "Esam") #Endotelial
FeaturePlot(GSE221226, features= "Pecam1") #Endotelial
FeaturePlot(GSE221226, features= "Lyz2") #Macrofago
FeaturePlot(GSE221226, features= "Cd68") #Macrofago
FeaturePlot(GSE221226, features= "Cnmd") #Inmune
FeaturePlot(GSE221226, features= "Cyp11a1") #lEY
FeaturePlot(GSE221226, features= "Igf1") #Estromal
FeaturePlot(GSE221226, features= "Acta2") #peritubular
GSE221226 <-  FindNeighbors(GSE221226, dims = 1:5)
GSE221226 <- FindClusters(GSE221226, resolution = c(0.3)) #11 COMUNIDADES
DimPlot(GSE221226, reduction = "umap", label = TRUE) #RESOLUTION.0.2 == CLUSTER 10
VlnPlot(GSE221226, features = "Ddx4")
VlnPlot(GSE221226, features = "Dazl")
cl_markers_GSE221226 <- FindAllMarkers(GSE221226,
                                       features = VariableFeatures(GSE221226))
cl_markers_GSE221226 <- cl_markers_GSE221226[cl_markers_GSE221226$avg_log2FC >= 1.5 & cl_markers_GSE221226$p_val_adj <= 0.05,]
df_GSE221226 <- GSE221226@meta.data
df = data.frame(RNA_snn_res.0.3 = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"), 
                cell_type = c("Sertoli cells", "Stromal cells", "Sertoli cells", "Sertoli cells", "Stromal cells", "Myoid cells", "Stromal cells", "Myoid cells",  "Germ cells", "Leydig", "Endothelial cells", "Endothelial cells"))
df_GSE221226$cells <- rownames(df_GSE221226)
df_GSE221226 <- data.frame (merge (df_GSE221226, df, by= "RNA_snn_res.0.3"))
GSE221226_UMAP <- data.frame (Embeddings(GSE221226, reduction = "umap"))
GSE221226_UMAP$cells <- rownames(GSE221226_UMAP)

df_GSE221226 <- merge(df_GSE221226, GSE221226_UMAP, by= "cells")
write.table(df_GSE221226, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSE221226")
write.table(cl_markers_GSE221226, file= "~/Figuras_Paper 2026/UMAP y Clustering/GSE221226_markers")

GSE221226_germinales <- GSE221226
dim (GSE221226_germinales)

samples <- list(                #Hay que volver a correrlo porque ahora los objetos seurat tienen los filtros aplicados
  GSM3744442_germinales =  GSM3744442_germinales,
  GSM3744443_germinales =  GSM3744443_germinales,   
  GSM3744444_germinales =  GSM3744444_germinales,
  GSM4185522_germinales = GSM4185522_germinales,
  GSM4185523_germinales =  GSM4185523_germinales, 
  GSM5708025_germinales = GSM5708025_germinales,
  GSM6468366_germinales = GSM6468366_germinales,
  GSM7135562_germinales = GSM7135562_germinales,
  SRR6129051_germinales = SRR6129051_germinales,
  SRR6129050_germinales = SRR6129050_germinales,
  GSE221226_germinales = GSE221226_germinales)

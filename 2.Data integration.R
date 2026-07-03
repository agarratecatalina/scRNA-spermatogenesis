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



#######Integración de datos 1########### Este va!!!!
#1) Hacemos merge de los objetos seurat

seurat.combined <- merge(
  GSM3744441_germinales,
  y = samples,
  add.cell.ids = c("GSM3744441_germinales",
                   "GSM3744442_germinales",
                   "GSM3744443_germinales",
                   "GSM3744444_germinales",
                   "GSM4185522_germinales",
                   "GSM4185523_germinales",
                   "GSM5708025_germinales",
                   "GSM6468366_germinales",
                   "GSM7135562_germinales",
                   "SRR6129051_germinales",
                   "SRR6129050_germinales",
                   "GSE221226_germinales"))

table (seurat.combined$orig.ident)

3)#control de calidad
df_seurat_combinated <- data.frame (seurat.combined@meta.data)
library (ggplot2)
ggplot(df_seurat_combinated, aes (x= df_seurat_combinated$orig.ident, y= df_seurat_combinated$nFeature_RNA)) + 
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.5, color = "darkblue")
ggplot(df_seurat_combinated, aes (x= df_seurat_combinated$nCount_RNA, y= df_seurat_combinated$nFeature_RNA, color= nFeature_RNA)) + 
  geom_point() 

dim (seurat.combined) #29720 12797
seurat.combined <- subset(seurat.combined, subset = nCount_RNA < 30000)
dim (seurat.combined) #29720 12700
DefaultAssay(seurat.combined) <- "RNA"
seurat.combined <- subset(seurat.combined, subset = nFeature_RNA > 20)

seurat.combined2 <- subset(
  seurat.combined,
  subset = (Sox9 == 0 )) #& 
              #Amh == 0 & 
              #Star == 0)) #SERTOLI
dim (seurat.combined2)
table (seurat.combined2$orig.ident)


seurat.combined2 <- subset(
  seurat.combined2,
  subset = (Cyp11a1 == 0)) #leydig

seurat.combined2 <- subset(
  seurat.combined2,
  subset = (Esam == 0 |
              Pecam1 == 0)) #ENDOTELIALES

seurat.combined2 <- subset(
  seurat.combined2,
  subset = (Acta2 == 0 )) #PERITUBULARES

seurat.combined2 <- subset(
  seurat.combined2,
  subset = (Lyz2 == 0 )) #Macrogagos

seurat.combined2 <- subset(
  seurat.combined2,
  subset = (Igf1 == 0)) #ESTROMALES


seurat.combined2 <- subset(
  seurat.combined2,
  subset = (Cnmd == 0)) #INMUNES

seurat.combined <- seurat.combined2
dim (seurat.combined)
meta <- seurat.combined@meta.data
meta$outlier <- meta$nCount_RNA > ave(
  meta$nCount_RNA,
  as.factor(meta$seurat_clusters),
  FUN = function(x) quantile(x, 0.95))
table (meta$outlier)
cells_keep <- rownames(meta)[!meta$outlier]
seurat.combined <- subset(seurat.combined, cells = cells_keep)

#2) Harmony
seurat.combined  <- NormalizeData(seurat.combined)
seurat.combined <- FindVariableFeatures(seurat.combined, nfeatures = 2000)
seurat.combined  <- ScaleData(seurat.combined)
seurat.combined  <-  RunPCA(seurat.combined, npcs = 100)
DimPlot(seurat.combined, reduction = "pca", group.by = "orig.ident") 
ElbowPlot(seurat.combined)


germinales_1000 <- IntegrateLayers(
  object = seurat.combined, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE)
DimPlot(germinales_1000, group.by = "orig.ident", reduction = "harmony")
ElbowPlot(germinales_1000)
dim (germinales_1000)

germinales <- RunUMAP(germinales, reduction = "harmony", dims = 1:5, reduction.name = "umap.harmony_pc5")
DimPlot(germinales, reduction = "umap.harmony_pc5",group.by = "orig.ident")
DimPlot(germinales, reduction = "umap.harmony_pc5", split.by = "orig.ident")

germinales <- RunUMAP(germinales, reduction = "harmony", dims = 1:8, reduction.name = "umap.harmony_pc8")
DimPlot(germinales, reduction = "umap.harmony_pc8",group.by = "orig.ident")

germinales <- RunUMAP(germinales, reduction = "harmony", dims = 1:10, reduction.name = "umap.harmony_pc10")
DimPlot(germinales, reduction = "umap.harmony_pc10",group.by = "orig.ident")
FeaturePlot(germinales, reduction = "umap.harmony_pc10", features = "nFeature_RNA")

germinales <- RunUMAP(germinales, reduction = "harmony", dims = 1:12, reduction.name = "umap.harmony_pc12")
DimPlot(germinales, reduction = "umap.harmony_pc12",group.by = "orig.ident")
FeaturePlot(germinales, reduction = "umap.harmony_pc12", features = "nFeature_RNA")

germinales <- RunUMAP(germinales, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony_pc15")
DimPlot(germinales, reduction = "umap.harmony_pc15",group.by = "orig.ident")
FeaturePlot(germinales, reduction = "umap.harmony_pc15", features = "nFeature_RNA")

germinales <- RunUMAP(germinales, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony_pc20")
DimPlot(germinales, reduction = "umap.harmony_pc20", group.by = "orig.ident")
FeaturePlot(germinales, reduction = "umap.harmony_pc20", features = "nFeature_RNA")

egerminales_1000 <- RunUMAP(germinales_1000, reduction = "harmony", dims = 1:25, reduction.name = "umap.harmony_pc25")
DimPlot(germinales_1000, reduction = "umap.harmony_pc25", group.by = "orig.ident")
FeaturePlot(germinales_1000, reduction = "umap.harmony_pc25", features = "nFeature_RNA")

germinales <- RunUMAP(germinales, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony_pc30")
DimPlot(germinales, reduction = "umap.harmony_pc30", group.by = "orig.ident")
FeaturePlot(germinales, reduction = "umap.harmony_pc30", features = "nFeature_RNA")

germinales <- RunUMAP(germinales, reduction = "harmony", dims = 1:35, reduction.name = "umap.harmony_pc35")
DimPlot(germinales, reduction = "umap.harmony_pc35", group.by = "orig.ident")
FeaturePlot(germinales, reduction = "umap.harmony_pc35", features = "nFeature_RNA")

germinales <- RunUMAP(germinales, reduction = "harmony", dims = 1:40, reduction.name = "umap.harmony_pc40")
DimPlot(germinales, reduction = "umap.harmony_pc40", group.by = "orig.ident")
FeaturePlot(germinales, reduction = "umap.harmony_pc40", features = "nFeature_RNA")

germinales <- RunUMAP(germinales, reduction = "harmony", dims = 1:50, reduction.name = "umap.harmony_pc50")
DimPlot(germinales, reduction = "umap.harmony_pc50", group.by = "orig.ident")
FeaturePlot(germinales, reduction = "umap.harmony_pc50", features = "nFeature_RNA")

FeaturePlot(germinales, features = "Nanos3",  reduction= "umap.harmony_pc20")
FeaturePlot(germinales, features = "Stra8",  reduction= "umap.harmony_pc15")
FeaturePlot(germinales, features = "Meiob",  reduction= "umap.harmony_pc15")
FeaturePlot(germinales, features = "Sycp3",  reduction= "umap.harmony_pc15")
FeaturePlot(germinales, features = "Acrv1",  reduction= "umap.harmony_pc20")
FeaturePlot(germinales, features = "Cers3",  reduction= "umap.harmony_pc15")
FeaturePlot(germinales, features = "Elovl2",  reduction= "umap.harmony_pc15")
FeaturePlot(germinales, features = "Tnp1",  reduction= "umap.harmony_pc20")

germinales <- FindNeighbors (germinales, reduction = "harmony", dims = 1:20)
germinales <- FindClusters(germinales, resolution = 0.5, cluster.name = "harmony_clusters")
germinales <- FindClusters(germinales, resolution = 0.2, cluster.name = "harmony_clusters")
germinales <- FindClusters(germinales, resolution = 1, cluster.name = "harmony_clusters")
DimPlot(germinales, reduction = "umap.harmony_pc20", label = TRUE)
DimPlot(germinales, reduction = "umap.harmony_pc20", label = TRUE, split.by = "orig.ident")
table (germinales$orig.ident, germinales$harmony_clusters)
DefaultAssay(germinales) <- "RNA" 
germinales <- JoinLayers(germinales, assay = "RNA")
cl_markers_germinales <- FindAllMarkers(germinales, 
                                        features = VariableFeatures(germinales))
cl_markers_germinales <- cl_markers_germinales[cl_markers_germinales$avg_log2FC >= 1.5 & cl_markers_germinales$p_val_adj <= 0.05 ,]
cl_markers_germinales2<- cl_markers_germinales[cl_markers_germin,etaales$avg_log2FC >= 2 & cl_markers_germinales$p_val_adj <= 0.05 ,]
cat (cl_markers_germinales[cl_markers_germinales$cluster == 8,7], sep= "\n")

germinales_finales <- subset (germinales, ident= c(0,1,2,3,4,5,6,7,8,10,11,12))

germinales_finales <- RunUMAP(germinales_finales, reduction = "harmony", dims = 1:10, reduction.name = "umap.harmony_pc10")





new_cluster_names <- c("Sertoli", "Estromales", "Sertoli", "Mioides", "Sertoli", "Spg indeferenciada", "Spg diferenciada", "Esp primario tardio", "Esp primario temprano", "Estromales", "Esp redonda", "Leydig", "Esp secundario", "Mioides", "Esp condensada", "Esp redonda", "Sertoli", "Endotelial", "Endotelial", "Macrofagos", "Endotelial")
names (new_cluster_names) <- levels (obj)
obj <- RenameIdents(obj, new_cluster_names)
Idents (obj)
table(Idents(obj), obj$orig.ident) #Veo cuantas celulas tengo en cada orig.
df_counts <- as.data.frame(table(cluster = Idents(obj), sample  = obj$orig.ident))
ggplot(df_counts, aes(x = sample, y = Freq, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  ylab("Porcentaje de células") +
  xlab("Muestra")

obj_subset <- subset(x = obj, idents = c("Spg indeferenciada", "Spg diferenciada", "Esp primario temprano", "Esp primario tardio", "Esp secundario", "Esp redonda", "Esp condensada"))
nuevo_orden <- c("Spg indeferenciada", "Spg diferenciada", "Esp primario temprano", "Esp primario tardio", "Esp secundario", "Esp condensada")


table(Idents(obj_subset), obj_subset$orig.ident) #Veo cuantas celulas tengo en cada orig.
df_counts <- as.data.frame(table(cluster = Idents(obj_subset), sample  = obj_subset$orig.ident))
ggplot(df_counts, aes(x = sample, y = Freq, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  ylab("Porcentaje de células") +
  xlab("Muestra")
DimPlot(obj_subset, reduction = "umap.harmony")
saveRDS(obj_subset, file = "C:/Users/Usuario/Desktop/Beca PosDoc/scRNAseq/obj_subset")

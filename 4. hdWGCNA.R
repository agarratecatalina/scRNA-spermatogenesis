# Install Bioconductor dependencies
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("WGCNA", "UCell", "GenomicRanges", "GeneOverlap"))
# Install hdWGCNA from GitHub
devtools::install_github("smorabit/hdWGCNA")


library (hdWGCNA)

#1. Genes de interes
features_2000 <- VariableFeatures(germinales) #2000 genes variables
germinales <- SetupForWGCNA (germinales,
               gene_select = "custom",
               gene_list = features_2000,
               wgcna_name = "scWGCNA_2000features")
GetWGCNAGenes (germinales) #Chequeamos si la lista de genes de interes se cargo bien

#2. Organizacion de las celulas en conjuntos o metacelulas
emb <- Embeddings(germinales, "umap.harmony_pc20")
rownames(emb) <- colnames(germinales)
germinales[["umap.harmony_pc20"]]@cell.embeddings <- emb

germinales <- MetacellsByGroups (germinales,
                                 group.by =  "harmony_clusters", #Construyo las metacelulas a partir de celulas de un mismo cluster
                                 reduction = "harmony", #Utilizo como partida las coordenadas del umap
                                 k = 20, #Por cada gen tendremos el promedio de 20 celulas
                                 max_shared = 5,
                                 ident.group = "harmony_clusters") 

#3. Con las metacelulas ya organizadas, debemos normalizar nuevamente
germinales <- NormalizeMetacells (germinales)

#4. Construimos una matriz de expresion de genes x metacelulas
germinales <- SetDatExpr(germinales,
                          assay = "RNA",
                          slot = "data")
#5. Seleccion del sodt treshold power
germinales <- TestSoftPowers (germinales)

#6. Construccion de la red de cooexpresion
germinales <- ConstructNetwork (germinales,
                                soft_power = 8,
                                  wgcna_name = "scWGCNA_2000features",
                                  overwrite_tom= TRUE,
                                tom_name = "INH")
PlotDendrogram(germinales)
module_palette <- c(
  grey      = "grey",
  brown     = "brown",
  red       = "red3",
  turquoise = "turquoise",
  blue      = "blue",
  pink      = "pink",
  green     = "green",
  black     = "black",
  yellow    = "yellow")
module_colors <- germinales@misc$scWGCNA_2000features$wgcna_net$colors
module_colors_custom <- module_palette[module_colors]
dend <- germinales@misc$scWGCNA_2000features$wgcna_net$dendrograms[[1]]

windowsFonts(Calibri = windowsFont("Calibri"))
par(family = "Calibri", cex = 1)
plotDendroAndColors(
       dend,
       module_colors_custom,
       groupLabels = "Module colors",
       dendroLabels = FALSE,
       hang = 0.03,
       addGuide = TRUE,
       guideHang = 0.05,
       main = "",
       cex.axis = 15 / 12,
       cex.lab = 15 / 12,
       cex.main = 15 / 12)

GetModules(germinales)
modulos <- GetModules(germinales)
blue_genes <- modulos [modulos$module == "blue",1] #fuerte en spg y va disminuyendo en espermatocitos
yellow_genes <- modulos [modulos$module == "yellow",1] #elogadas
black_genes <- modulos [modulos$module == "black",1] #redondas
green_genes <- modulos [modulos$module == "green",1] #spII y redondas
red_genes <- modulos [modulos$module == "red",1] #SpI
pink_genes <- modulos [modulos$module == "pink",1] #entre spI-spII
brown_genes <- modulos [modulos$module == "brown",1] #spg inicial
turquoise_genes <- modulos [modulos$module == "turquoise",1]  #spg y spI

#7.
germinales <- ModuleEigengenes(germinales)
MEs <- GetMEs(germinales)
head(MEs)
germinales@meta.data <- cbind(germinales@meta.data, MEs)
MEs <- MEs[, c(1:7,9)]

ME_cor <- cor(MEs)
ME_diss <- 1 - ME_cor
hc <- hclust(as.dist(ME_diss), method = "average")
plot(hc)

cor_ME <- cor(MEs)
library(pheatmap)
rownames(cor_ME) <- c("Yellow", "Black", "Green", "Red", "Pink", "Brown", "Turquoise", "Blue")
colnames(cor_ME) <- c("Yellow", "Black", "Green", "Red", "Pink", "Brown", "Turquoise", "Blue")
pheatmap(cor_ME,
         fontsize_row = 12,
         fontsize_col = 12,
         angle_col = 45)

#8. Idebtufucar genes hub
germinales <- ModuleConnectivity(germinales)
hub_df <- GetHubGenes(germinales, 50)

#9. TOM
matriz_TOM <- GetTOM(germinales)



###############################
library (clusterProfiler)
library (org.Mm.eg.db)

genes_symbol_lista <- list (
  blue = blue_genes,
  brown = brown_genes,
  yellow = yellow_genes,
  red = red_genes,
  pink = pink_genes,
  green = green_genes,
  black = black_genes,
  turquoise = turquoise_genes
)

genes_entrez_lista <- lapply(genes_symbol_lista, function(genes) {
  df_conversion <- bitr(
    genes, 
    fromType = "SYMBOL", 
    toType   = "ENTREZID", 
    OrgDb    = org.Mm.eg.db, 
    drop     = TRUE
  )
  return(df_conversion$ENTREZID)
})

comparacion_GO_MF <- compareCluster(
  geneClusters  = genes_entrez_lista,
  fun           = "enrichGO",
  OrgDb         = org.Mm.eg.db,
  ont           = "MF", # BP = Proceso Biológico
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# 4. Visualizar el resultado
dotplot(comparacion_GO_MF, title = "Comparación GO (Entrez ID)")


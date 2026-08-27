library(SeuratWrappers)
library (monocle3)
cds <- as.cell_data_set(obj_subset) #transformo mi objeto seurat en un objeto monocle
# Clusterizar células (necesario para pseudotime)
cds <- cluster_cells(cds)
# Aprender el grafo de trayectoria
cds <- learn_graph(cds)
# Ordenar células a lo largo del pseudotime
cds <- order_cells(cds)
# Graficar pseudotime
plot_cells(cds, color_cells_by = "pseudotime")
# Evaluamos dentro de nuestros genes lipidicos cuales presentan variacion pseudotemporal
cds_subset <- cds[genes_presentes, ]
test_res <- graph_test(cds_subset, 
                    neighbor_graph = "principal_graph")
# Ordenar el data frame generado de acuerdo al q_value
test_res <- test_res[order(test_res$q_value), ]
# Filtramos los genes que son significativos
sig_genes <- subset(test_res, q_value < 0.05)
genes <- rownames(sig_genes) #Nos quedamos con los nombres de los genes
# matriz genes x células
expr_mat <- exprs(cds_subset[genes, ]) #Volvemos a hacer un subset para quedarnos solo con los genes que fueron significativos
# ordenar células por pseudotiempo
pseudo <- pseudotime(cds_subset) #Guardamos la información pseudotempral en un vector.
cells_order <- order(pseudo)  #VER BIEN ESTO
expr_mat <- expr_mat[, cells_order] #VER BIEN ESTO


# opcional: log-transformar
expr_mat <- log1p(expr_mat)
#Escalar datos:
expr_mat_scaled <- t(scale(t(expr_mat)))
#K-means:
k <- 4  # número de clusters que quieras
clusters <- kmeans(expr_mat_scaled, centers = k)$cluster
table(clusters)
#H-clust:
dist_mat <- dist(expr_mat_scaled)
hc <- hclust(dist_mat)
plot(hc)


# Matriz de expresión (genes x células)
expr_mat <- counts(cds_subset)
# Fracción de células que expresan cada gen
gene_fraction <- Matrix::rowSums(expr_mat > 0) / ncol(expr_mat)
# Seleccionar genes expresados en >=25% de células
genes_presentes2 <- names(gene_fraction[gene_fraction >= 0.25])
cds_subset <- cds[genes_presentes, ]
test_res <- graph_test(cds_subset, 
                       neighbor_graph = "principal_graph", 
                       cores = 4)
# Ordenar por q_value
test_res <- test_res[order(test_res$q_value), ]
# Filtrar significativos
sig_genes <- subset(test_res, q_value < 0.05)

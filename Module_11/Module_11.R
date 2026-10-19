#########################
### trajectory analysis
## Pseudotime example ##
#########################
BLCA_M11 <- readRDS("BLCA_M11.RDS")

library(Seurat)
library(VennDiagram)
library(tidyverse)
library(ggpubr)
library(plotly)
library(cowplot)
library(SingleCellExperiment)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(SeuratWrappers)
library(monocle3)
######################################
######################################
## https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02341-y
remotes::install_github('cole-trapnell-lab/monocle3', force = TRUE)

# Convert to monocle3 cell_data_set
cells.cds <- as.cell_data_set(BLCA_M11)

# # to get cell meta.data ....
# colData(cells.cds)

# # to get gene meta.data
# fData(cells.cds)
# rownames(fData(cells.cds))[1:10]
# fData(cells.cds)$gene_short_name <- rownames(fData(cells.cds))

# # to get counts
# counts(cells.cds)

# Re-cluster (required for trajectories in monocle3)
cells.cds <- cluster_cells(cells.cds)

# Learn the trajectory
cells.cds <- learn_graph(cells.cds)

# Plot the trajectory
####################################################################
names(BLCA_M11@assays) <- c("SCT","RNA")
BLCA_M11@active.assay = 'RNA'

rowData(cells.cds)$gene_name <- rownames(cells.cds)
rowData(cells.cds)$gene_short_name <- rowData(cells.cds)$gene_name
###################################################################
rami1 <- plot_cells(cells.cds, color_cells_by = "seurat_clusters", 
                    show_trajectory_graph = TRUE)
ggsave(filename = "QC181.png", plot = rami1, width = 8, height = 6)
###################################
###################################
plot_cells(cells.cds, 
           color_cells_by = "cell_type", 
           show_trajectory_graph = TRUE, genes = "TIMP1")

rami2 <- plot_cells(cells.cds, 
                    color_cells_by = "cell_type",
                    show_trajectory_graph = TRUE,  genes = "TIMP1")
ggsave(filename = "QC182.png", plot = rami2, width = 8, height = 6)
##################################
##################################
rami1 <- plot_cells(cells.cds, color_cells_by = "seurat_clusters",
                    reduction_method = "UMAP", 
                    label_cell_groups = FALSE,
                    label_leaves = FALSE,
                    show_trajectory_graph = FALSE, 
                    genes = "TIMP1")
ggsave(filename = "QC183.png", plot = rami1, width = 8, height = 6)
##################################
##################################
features=c("TIMP1", "CD63", "MMP9", 
           "MMP14", "IL6", "FN1")

rami1 <- plot_cells(cells.cds, color_cells_by = "seurat_clusters",
                    reduction_method = "UMAP", 
                    label_cell_groups = FALSE,
                    label_leaves = FALSE,
                    show_trajectory_graph = FALSE, 
                    genes = features)
ggsave(filename = "QC184.png", plot = rami1, width = 12, height = 8)
###########################################
###########################################
# Find the pseudotime order
DimPlot(BLCA_M11, group.by = "Phase")
cells.cds <- order_cells(cells.cds, reduction_method = "UMAP")

# Plot pseudotime
rami5 <- plot_cells(cells.cds, color_cells_by = "pseudotime")
rami6 <- plot_cells(cells.cds, color_cells_by = "pseudotime", genes = features)
rami7 <- plot_cells(cells.cds, color_cells_by = "pseudotime", genes = "TIMP1")

ggsave(filename = "QC185.png", plot = rami5, width = 8, height = 6)
ggsave(filename = "QC186.png", plot = rami6, width = 12, height = 8)
ggsave(filename = "QC187.png", plot = rami7, width = 8, height = 6)
#######################################################################
#######################################################################
# Add pseudotime score back into the dataset
BLCA_M11$pseudotime_score <- cells.cds@principal_graph_aux@listData$UMAP$pseudotime
BLCA_M11$pseudotime_score[BLCA_M11$pseudotime_score == Inf] <- NA

# FeaturePlot with Pseudotime
FeaturePlot(BLCA_M11, features = c("pseudotime_score"))
FeaturePlot(BLCA_M11, features = c("pseudotime_score")) + DimPlot(BLCA_M11, group.by = "Phase")

DimPlot(BLCA_M11, group.by = "Phase") + FeaturePlot(BLCA_M11, features = c("pseudotime_score"))

rami1 <- DimPlot(BLCA_M11, group.by = "Phase")
rami2 <- FeaturePlot(BLCA_M11, features = c("pseudotime_score"))
rami3 <- (rami1 + rami2) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC188.png", width = 12, height = 6)

rami1 <- FeaturePlot(BLCA_M11, features = "TIMP1") 
rami2 <- FeaturePlot(BLCA_M11, features = c("pseudotime_score"))

rami3 <- ((rami1) | (rami2)) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC189.png", plot = rami3, width = 12, height = 6)
#######################################
#######################################
rami1 <- FeaturePlot(BLCA_M11, features = "TIMP1", reduction = "tsne") 
rami2 <- FeaturePlot(BLCA_M11, features = c("pseudotime_score"), reduction = "tsne")

rami3 <- ((rami1) | (rami2))
ggsave(filename = "QC190.png", plot = rami3, width = 12, height = 6)
###################################################################
saveRDS(BLCA_M11, file = "BLCA_M12.RDS")
BLCA_M12 <- readRDS("BLCA_M12.RDS")
#######################################
#######################################
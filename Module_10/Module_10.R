## Load libraries ##
library(Seurat)
library(VennDiagram)
library(tidyverse)
library(ggpubr)
library(plotly)
library(cowplot)
library(SingleCellExperiment)
library(patchwork)
library(pheatmap)
library(enrichR)
library(RColorBrewer)
library(SingleR)
##############################
##############################
### Read your data
BLCA_M11 <- readRDS("BLCA_M11.RDS")
#########################################
# Helper function for plotting enrichr results
plot_eres <- function(eres_name, eres_list, n = 10) {
  eres_list[[eres_name]] %>%
    top_n(n = n, wt = -log10(Adjusted.P.value)) %>%
    arrange(-log10(Adjusted.P.value)) %>%
    mutate(Term = factor(Term, levels = Term)) %>%
    ggplot(mapping = aes(x = Term, y = -log10(Adjusted.P.value), 
                         fill = Combined.Score)) +
    geom_bar(stat = "identity") +
    ggpubr::rotate() +
    theme_bw(base_size = 16) +
    rremove("ylab") +
    labs(title = eres_name)}
########################################
########################################
# Define color palettes
WtOrRd_pal <- colorRampPalette(c("lightblue", brewer.pal(9, "OrRd")))(100)
# Get enrichr databases
dbs <- listEnrichrDbs()
# Set random seed
set.seed(42) 
########################################
# Dim plot of UMAP
Idents(BLCA_M11) <- "seurat_clusters"

rami1 <- DimPlot(BLCA_M11, group.by = "seurat_clusters", 
                 label = TRUE, label.size = 6)
rami2 <- FeaturePlot(BLCA_M11, features = "TIMP1", 
                     label = TRUE, label.size = 6) ## 
FeaturePlot(BLCA_M11, features = "TIMP1", label = TRUE) ## 
rami3 <- (rami1 + rami2) & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC172.png", plot = rami3, width = 20, height = 10)
#########################################
## What's up with unknown cells?
DimPlot(BLCA_M11, group.by = "seurat_clusters") + 
  FeaturePlot(BLCA_M11, features = "nFeature_RNA")
p1 <- DimPlot(BLCA_M11, group.by = "seurat_clusters", reduction = "tsne", 
              label = TRUE, label.size = 4.5)
p2 <- FeaturePlot(BLCA_M11, features = "TIMP1", reduction = "tsne", 
                  label = TRUE, label.size = 4.5)
p3 <- p1 + p2
ggsave(filename = "QC173.png", plot = p3, width = 20, height = 10)
##########################################
Rami1 <- DimPlot(BLCA_M11, group.by = "seurat_clusters",
                 reduction = "umap", label = TRUE,
                 label.size = 6) + FeaturePlot(BLCA_M11, 
                                               features = "nFeature_RNA")
rami2 <- Rami1 & labs(x = "UMAP_1", y = "UMAP_2")
ggsave(filename = "QC174.png", plot = rami2, width = 20, height = 10)
######################################################
######################################################
## Find cluster markers ##
BLCA.markers <- FindAllMarkers(BLCA_M11, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.5)

sig_markers <- BLCA.markers %>% 
  filter(p_val_adj < 0.05) %>%
  write_csv(file = "BLCA_markers.csv")

top_markers <- sig_markers %>% 
  group_by(cluster) %>% 
  top_n(n = 50, wt = avg_log2FC) %>%
  write_csv(file = "top_markers.csv")
#########################################
#########################################
# Cluster 0 enrichr: https://maayanlab.cloud/Enrichr/
to_check <- c("Human_Gene_Atlas", "GO_Cellular_Component_2023",
              "GO_Molecular_Function_2023", "GO_Biological_Process_2023")

# Pull the Cluster 0 genes as a vector
cluster1_genes <- sig_markers %>%
  filter(cluster == 6) %>%
  pull(gene)

cluster2_genes <- sig_markers %>%
  filter(cluster == 11) %>%
  pull(gene)

cluster3_genes <- sig_markers %>%
  filter(cluster == 18) %>%
  pull(gene)

cluster4_genes <- sig_markers %>%
  filter(cluster == 19) %>%
  pull(gene)

cluster5_genes <- sig_markers %>%
  filter(cluster == 24) %>%
  pull(gene)
#####################################
#####################################
# Run through enrichr
cluster1_eresList <- enrichr(cluster1_genes, databases = to_check)
cluster2_eresList <- enrichr(cluster2_genes, databases = to_check)
cluster3_eresList <- enrichr(cluster3_genes, databases = to_check)
cluster4_eresList <- enrichr(cluster4_genes, databases = to_check)
cluster5_eresList <- enrichr(cluster5_genes, databases = to_check)

# Get the top hits from each and plot them
library(ggpubr)
?cowplot

plotList <- lapply(names(cluster1_eresList), plot_eres, eres_list = cluster1_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh", ncol = 2) +
  patchwork::plot_annotation(title = "Cluster [6] Enrichr Analysis", 
                             theme = theme(title = element_text(size = 20))) 
ggsave(filename = "QC175.png", width = 25, height = 15)
###############################
plotList <- lapply(names(cluster2_eresList), plot_eres, eres_list = cluster2_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh", ncol = 2) +
  patchwork::plot_annotation(title = "Cluster [11] Enrichr Analysis", 
                             theme = theme(title = element_text(size = 20))) 
ggsave(filename = "QC176.png", width = 25, height = 15)
################################
plotList <- lapply(names(cluster3_eresList), plot_eres, eres_list = cluster3_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh", ncol = 2) +
  patchwork::plot_annotation(title = "Cluster [18] Enrichr Analysis", 
                             theme = theme(title = element_text(size = 20))) 
ggsave(filename = "QC177.png", width = 25, height = 15)
#################################
plotList <- lapply(names(cluster4_eresList), plot_eres, eres_list = cluster4_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh", ncol = 2) +
  patchwork::plot_annotation(title = "Cluster [19] Enrichr Analysis", 
                             theme = theme(title = element_text(size = 20))) 
ggsave(filename = "QC178.png", width = 25, height = 15)
#################################
plotList <- lapply(names(cluster5_eresList), plot_eres, eres_list = cluster5_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh", ncol = 2) +
  patchwork::plot_annotation(title = "Cluster [24] Enrichr Analysis", 
                             theme = theme(title = element_text(size = 20))) 
ggsave(filename = "QC179.png", width = 25, height = 15)
####################################
####################################
### Batch ..... *******************

## Find markers of cluster 0 vs cluster 2 ##
cluster18vs19.markers <- FindMarkers((BLCA_M11),
                                    ident.1 = 18, 
                                    ident.2 = 19,
                                    logfc.threshold = .25,
                                    min.pct = 0.25)
# Pull the Cluster 8 vs 14 genes as a vector
cluster8_genes <- cluster18vs19.markers %>%
  filter(avg_log2FC > 0 & p_val_adj < 0.05) %>%
  rownames_to_column("gene") %>%
  pull(gene)
# Run through enrichr
enrichR::genes790

to_check <- c("Human_Gene_Atlas", "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
              "KEGG_2019_Human", "MSigDB_Hallmark_2020")

cluster0_eresList <- enrichr(cluster8_genes, databases = to_check)

# Get the top hits from each and plot them
plotList <- lapply(names(cluster0_eresList), plot_eres, eres_list = cluster0_eresList)
cowplot::plot_grid(plotlist = plotList, labels = "AUTO",
                   label_size = 18, align = "vh") +
  patchwork::plot_annotation(title = "Cluster 18 vs 19 Enrichr Analysis", 
                             theme = theme(title = element_text(size = 20))) 
ggsave(filename = "QC180.png", width = 25, height = 15)
#########################################
#########################################
install.packages("enrichR")
library(enrichR)

dbs <- listEnrichrDbs()
head(dbs)

dbs1 <- c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", 
          "GO_Biological_Process_2023", "KEGG_2021_Human",
          "GTEx_Tissues_V8_2023", "Cancer_Cell_Line_Encyclopedia")

enriched <- enrichr(genes = c("TIMP1", "MMP9"), databases = dbs1)
enriched <- enrichr(genes = "TIMP1", databases = dbs1)

printEnrich(enriched, "output.txt" , 
            #sep = "\t", 
            columns = c(1:9))
enriched[["GO_Biological_Process_2023"]]

plotEnrich(enriched[[1]], showTerms = 20, numChar = 40,
           y = "Count", orderBy = "P.value")
########################################################
########################################################
saveRDS(BLCA_M11, file = "BLCA_M11.RDS")
BLCA_M11 <- readRDS("BLCA_M11.RDS")
#########################################
#########################################
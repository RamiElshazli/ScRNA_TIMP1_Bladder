# TIMP1–BLCA | Reviewer 1, Comment 3 | Annotation validation (NO AnnQ)
# Input: BLCA_M7.RDS (Seurat object). R >= 4.3 recommended; Seurat >= 5.
# Outputs: CSV tables, PDF plots, and analysis log in reviewer1_comment3_outputs/
# This script does NOT infer malignancy, prove stemness, or perform CNV inference.

options(stringsAsFactors = FALSE)
INPUT_RDS <- "BLCA_M7.RDS"  # CHANGE to full path if needed
OUTDIR <- "reviewer1_comment3_outputs"
RUN_SINGLER_IF_MISSING <- TRUE
RUN_SECOND_REFERENCE <- FALSE  # optional reference sensitivity: BlueprintEncodeData
MIN_CELLS_FOR_MARKER_TEST <- 10L
set.seed(42)

required <- c("Seurat", "SeuratObject", "ggplot2", "dplyr", "tidyr", "Matrix", "SingleCellExperiment", "SingleR", "celldex", "scater", "patchwork")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Install missing packages before running: ", paste(missing, collapse=", "),
  "\nCRAN: install.packages(c('Seurat','SeuratObject','ggplot2','dplyr','tidyr','Matrix','patchwork'))",
  "\nBioconductor: BiocManager::install(c('SingleCellExperiment','SingleR','celldex','scater'))")

library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
dir.create(OUTDIR, showWarnings=FALSE, recursive=TRUE)
logfile <- file.path(OUTDIR, "analysis_log.txt")
logline <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n", file=logfile, append=TRUE)
writecsv <- function(x, name) utils::write.csv(x, file.path(OUTDIR, name), row.names=FALSE, na="")
if (!file.exists(INPUT_RDS)) stop("RDS not found: ", normalizePath(INPUT_RDS, mustWork=FALSE))
obj <- readRDS(INPUT_RDS)
if (!inherits(obj, "Seurat")) stop("Input must be a Seurat object")
logline("Loaded", INPUT_RDS, "cells", ncol(obj), "genes", nrow(obj))
logline("Metadata columns:", paste(colnames(obj@meta.data), collapse=", "))

# ---- Metadata identification: review mapping file and edit if auto-detection is wrong ----
md <- obj@meta.data
pick_col <- function(candidates, cols) {
  exact <- candidates[candidates %in% cols]
  if (length(exact)) return(exact[[1]])
  for (x in candidates) { z <- grep(x, cols, ignore.case=TRUE, value=TRUE); if(length(z)) return(z[[1]]) }
  NA_character_
}
SPECIMEN_COL <- pick_col(c("specimen", "sample_id", "sample", "orig.ident", "Sample", "patient_id"), names(md))
CELLTYPE_COL <- pick_col(c("cell_type", "celltype", "CellType", "cell_type_manual", "annotation", "cell_types", "predicted_celltype"), names(md))
SOURCE_COL <- pick_col(c("source_study", "study", "accession", "dataset", "source"), names(md))
DONOR_COL <- pick_col(c("donor_id", "donor", "patient_id", "patient"), names(md))
CLUSTER_COL <- pick_col(c("seurat_clusters", "cluster"), names(md))
writecsv(data.frame(field=c("specimen","celltype","source","donor","cluster"),
  detected_column=c(SPECIMEN_COL,CELLTYPE_COL,SOURCE_COL,DONOR_COL,CLUSTER_COL)), "00_metadata_column_mapping_CHECK.csv")
if (is.na(SPECIMEN_COL)) stop("Cannot identify specimen column. Inspect metadata and set SPECIMEN_COL manually.")
if (is.na(CELLTYPE_COL)) {
  logline("WARNING: No cell-type metadata detected; using Seurat Idents as provisional annotation")
  md$provisional_identity <- as.character(Idents(obj)); CELLTYPE_COL <- "provisional_identity"
}
if (is.na(CLUSTER_COL)) {md$cluster_from_idents <- as.character(Idents(obj)); CLUSTER_COL <- "cluster_from_idents"}
md$specimen_for_analysis <- as.character(md[[SPECIMEN_COL]])
md$celltype_for_analysis <- as.character(md[[CELLTYPE_COL]])
md$cluster_for_analysis <- as.character(md[[CLUSTER_COL]])
md$source_for_analysis <- if (!is.na(SOURCE_COL)) as.character(md[[SOURCE_COL]]) else NA_character_
md$donor_for_analysis <- if (!is.na(DONOR_COL)) as.character(md[[DONOR_COL]]) else NA_character_
# If source missing, use only the specimen IDs P1–P10, N1–N3, U1–U8 from this project's mapping.
if (all(is.na(md$source_for_analysis))) {
  s <- toupper(gsub("[^A-Za-z0-9]", "", md$specimen_for_analysis))
  study_map <- c(P1="10xGenomics",P2="GSE130001",P3="GSE130001",N1="GSE129845",N2="GSE129845",N3="GSE129845",C1="GSE129845",C2="GSE129845",C3="GSE129845")
  study_map <- c(study_map, setNames(rep("GSE267718",7),paste0("P",4:10)), setNames(rep("GSE267718",8),paste0("U",1:8)))
  md$source_for_analysis <- unname(study_map[s]);
  logline("Source study derived from project specimen mapping; verify sample IDs")
}
writecsv(data.frame(specimen=sort(unique(md$specimen_for_analysis)),
  n_cells=as.integer(table(md$specimen_for_analysis)[sort(unique(md$specimen_for_analysis))])),
  "01_specimen_inventory.csv")
writecsv(as.data.frame(md[,c("specimen_for_analysis","celltype_for_analysis","cluster_for_analysis","source_for_analysis","donor_for_analysis"),drop=FALSE]) %>%
  mutate(cell_barcode=rownames(md)), "02_cell_metadata.csv")

# ---- Original NON-IMPUTED RNA counts; Seurat v5 layers are joined if needed ----
if (!"RNA" %in% Assays(obj)) stop("RNA assay absent. Do not substitute MAGIC or integrated expression")
DefaultAssay(obj) <- "RNA"
obj <- tryCatch(JoinLayers(obj, assay="RNA"), error=function(e) {logline("JoinLayers warning:",conditionMessage(e));obj})
counts <- tryCatch(GetAssayData(obj, assay="RNA", layer="counts"), error=function(e) NULL)
if (is.null(counts)) stop("Unable to retrieve joined RNA counts. Inspect assay layers and join original RNA count layers manually.")
if (ncol(counts) != nrow(md)) stop("Count matrix and metadata dimensions do not match")
counts <- counts[,rownames(md),drop=FALSE]
logline("Raw RNA counts matrix:", nrow(counts), "genes x", ncol(counts), "cells")

# ---- Marker panels: include positive and potential exclusion/contamination markers ----
marker_panels <- list(
  Urothelial_epithelial=c("EPCAM","KRT8","KRT18","KRT19","UPK1A","UPK2","KRT5","KRT14","TP63"),
  Fibroblast=c("COL1A1","COL1A2","DCN","LUM","PDGFRA","FAP","ACTA2"),
  Endothelial=c("PECAM1","VWF","KDR","EMCN","CLDN5"),
  Myeloid=c("LST1","LYZ","CD68","CSF1R","C1QA","C1QB","S100A8","S100A9"),
  T_NK=c("CD3D","CD3E","TRAC","NKG7","GNLY","KLRD1"),
  B_Plasma=c("MS4A1","CD79A","CD79B","JCHAIN","MZB1"),
  Stem_associated_not_proof=c("KRT5","KRT14","KRT17","TP63","ITGA6","SOX2"),
  Neuronal_like=c("TUBB3","RBFOX3","SNAP25","SYT1","ELAVL4"),
  Myeloid_progenitor_like=c("CD34","KIT","GATA2","MPO","ELANE"),
  Exclusion_doublet_review=c("PTPRC","EPCAM","COL1A1","PECAM1","HBB","MKI67")
)
markers <- unique(unlist(marker_panels)); present <- intersect(markers, rownames(counts))
writecsv(data.frame(marker=markers, present_in_raw_RNA=markers %in% present,
  panels=vapply(markers, function(x) paste(names(Filter(function(v) x %in% v, marker_panels)), collapse="; "), character(1))),
  "03_marker_availability.csv")
if (!length(present)) stop("None of the expected human gene symbols are present; inspect gene IDs")

# For large datasets, compute per-celltype and per-specimen marker detection efficiently from sparse counts.
calculate_marker_stats <- function(groups, group_label) {
  ix <- split(seq_along(groups), as.character(groups)); ix <- ix[lengths(ix)>0]
  res <- lapply(names(ix), function(g) {
    ids <- ix[[g]]; sub <- counts[present,ids,drop=FALSE]
    data.frame(group=g, group_type=group_label, marker=present, n_cells=length(ids),
      n_positive=as.numeric(Matrix::rowSums(sub>0)),
      pct_positive=100*as.numeric(Matrix::rowSums(sub>0))/length(ids),
      mean_raw_umi=as.numeric(Matrix::rowMeans(sub)))
  }); dplyr::bind_rows(res)
}
celltype_stats <- calculate_marker_stats(md$celltype_for_analysis,"cell_type")
specimen_stats <- calculate_marker_stats(md$specimen_for_analysis,"specimen")
writecsv(celltype_stats,"04_marker_by_celltype.csv")
writecsv(specimen_stats,"05_marker_by_specimen.csv")
# Cell-type dot plot: size = fraction detected, color = mean raw UMI (not normalized).
plotdata <- celltype_stats %>% mutate(marker=factor(marker,levels=rev(present)))
p <- ggplot(plotdata,aes(x=group,y=marker,size=pct_positive,color=log1p(mean_raw_umi))) +
  geom_point() + scale_size(range=c(0.2,6),limits=c(0,100)) +
  labs(x="Provisional cell-type annotation",y="Marker gene",size="Detected (%)",color="log(1 + mean raw UMI)",
       title="Lineage-marker evidence (non-imputed RNA counts)") +
  theme_bw(base_size=10)+theme(axis.text.x=element_text(angle=55,hjust=1))
ggsave(file.path(OUTDIR,"Figure_S4_marker_dotplot.pdf"),p,width=max(11,0.6*length(unique(plotdata$group))+5),height=12,limitsize=FALSE)

# ---- Cell-type counts by specimen / source / donor (if donor metadata exists) ----
ct <- md %>% count(specimen_for_analysis,source_for_analysis,donor_for_analysis,celltype_for_analysis,name="n_cells") %>%
  group_by(specimen_for_analysis) %>% mutate(within_specimen_pct=100*n_cells/sum(n_cells)) %>% ungroup()
writecsv(ct,"Table_S5_celltype_counts_by_specimen.csv")
heat <- ggplot(ct,aes(x=specimen_for_analysis,y=celltype_for_analysis,fill=within_specimen_pct)) +
  geom_tile(color="white",linewidth=.1)+scale_fill_gradient(low="white",high="#20609b")+
  labs(x="Specimen",y="Provisional cell type",fill="Within specimen (%)",
       title="Cell-type composition by specimen (descriptive only)") +
  theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=60,hjust=1))
ggsave(file.path(OUTDIR,"Figure_S5_celltype_by_specimen.pdf"),heat,width=12,height=7)

# ---- SingleR confidence: rerun if original prediction scores unavailable ----
# Rerun on non-imputed RNA, not integrated / MAGIC assay.
# WARNING: Human Primary Cell Atlas is not bladder-specific and may mislabel unusual states.
run_singler <- function(ref_object, ref_name) {
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts))
  # SingleR uses test logcounts if available; library-normalize original raw RNA counts.
  sce <- scater::logNormCounts(sce)
  ref_labels <- as.character(SummarizedExperiment::colData(ref_object)$label.main)
  if (!"logcounts" %in% SummarizedExperiment::assayNames(ref_object))
    stop("Reference lacks logcounts")
  pred <- SingleR::SingleR(test=sce, ref=ref_object, labels=ref_labels,
    assay.type.test="logcounts", assay.type.ref="logcounts", de.method="classic")
  scores <- as.matrix(pred$scores)
  # Score margin is a DESCRIPTIVE metric; it is NOT a calibrated posterior probability.
  top2 <- t(apply(scores,1,function(z) sort(z,decreasing=TRUE)[seq_len(min(2,length(z)))]))
  margin <- if(ncol(top2)>=2) top2[,1]-top2[,2] else rep(NA_real_,nrow(top2))
  out <- data.frame(cell_barcode=rownames(pred),reference=ref_name,
    predicted_label=as.character(pred$labels),
    pruned_label=as.character(pred$pruned.labels),
    delta_next=as.numeric(pred$delta.next),
    top_score=as.numeric(top2[,1]),score_margin_top2=as.numeric(margin))
  out <- cbind(out,md[out$cell_barcode,c("specimen_for_analysis","source_for_analysis","donor_for_analysis","celltype_for_analysis"),drop=FALSE])
  writecsv(out,paste0("06_SingleR_confidence_",ref_name,".csv"))
  summary <- out %>% group_by(celltype_for_analysis,specimen_for_analysis,source_for_analysis) %>%
    summarise(n=n(),pct_pruned=100*mean(is.na(pruned_label)),
      median_delta_next=median(delta_next,na.rm=TRUE),
      median_score_margin=median(score_margin_top2,na.rm=TRUE),.groups="drop")
  writecsv(summary,paste0("07_SingleR_confidence_summary_",ref_name,".csv"))
  fig <- ggplot(out,aes(x=celltype_for_analysis,y=delta_next))+
    geom_boxplot(outlier.shape=NA)+coord_flip()+
    labs(title=paste("SingleR annotation confidence:",ref_name),x="Provisional cell type",y="SingleR delta.next (not a probability)")+
    theme_bw(base_size=10)
  ggsave(file.path(OUTDIR,paste0("Figure_S6_SingleR_confidence_",ref_name,".pdf")),fig,width=10,height=7)
  out
}
if (RUN_SINGLER_IF_MISSING) {
  logline("Starting SingleR with HumanPrimaryCellAtlasData; may take substantial RAM/time")
  hpca <- celldex::HumanPrimaryCellAtlasData()
  hpca_pred <- run_singler(hpca,"HPCA")
  if (RUN_SECOND_REFERENCE) {
    logline("Starting optional alternative reference: BlueprintEncodeData")
    blue <- celldex::BlueprintEncodeData()
    blue_pred <- run_singler(blue,"BlueprintEncode")
    comp <- hpca_pred %>% select(cell_barcode,HPCA_label=predicted_label) %>%
      inner_join(blue_pred %>% select(cell_barcode,Blueprint_label=predicted_label),by="cell_barcode")
    writecsv(comp,"08_reference_comparison_raw_labels.csv")
    logline("CAUTION: Different reference ontologies; raw label disagreement is not directly equivalent to annotation error")
  }
}

# ---- Rare / unusual annotation reassessment (descriptive, not proof of lineage) ----
rare_pattern <- "stem|neuron|neural|progenitor|myeloid progenitor|cmp"
rare_types <- unique(md$celltype_for_analysis[grepl(rare_pattern,md$celltype_for_analysis,ignore.case=TRUE)])
if(length(rare_types)) {
  rare_ct <- ct %>% filter(celltype_for_analysis %in% rare_types)
  writecsv(rare_ct,"09_rare_celltype_specimen_support.csv")
  rare_marker <- celltype_stats %>% filter(group %in% rare_types)
  writecsv(rare_marker,"10_rare_celltype_marker_support.csv")
  # Marker co-detection can indicate unusual mixed signatures; not definitive doublets.
  rare_ids <- which(md$celltype_for_analysis %in% rare_types)
  candidates <- intersect(c("EPCAM","PTPRC","COL1A1","PECAM1","TUBB3","RBFOX3","CD34","KIT","MPO","KRT5","KRT14"),rownames(counts))
  if(length(candidates)) {
    x <- as.matrix(counts[candidates,rare_ids,drop=FALSE]>0)
    support <- data.frame(cell_barcode=colnames(counts)[rare_ids],
      proposed_cell_type=md$celltype_for_analysis[rare_ids],specimen=md$specimen_for_analysis[rare_ids],
      source_study=md$source_for_analysis[rare_ids],x,check.names=FALSE)
    writecsv(support,"11_rare_cell_marker_codetection.csv")
  }
} else logline("No rare-cell annotations matched patterns; inspect labels in metadata")

# ---- Final audit ----
logline("Completed. Output directory:",normalizePath(OUTDIR))
logline("WARNING: No AnnQ, CNV inference, pathology, doublet calling, or formal reference coverage benchmark performed.")
logline("WARNING: SingleR score margins and delta.next are not posterior probabilities.")
logline("WARNING: If donor ID is missing, specimen-level counts MUST NOT be called per-donor counts.")
cat("Finished. See ",normalizePath(OUTDIR),"\n")

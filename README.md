# TIMP1 expression heterogeneity in bladder carcinoma

**Repository:** [ScRNA_TIMP1_Bladder](https://github.com/RamiElshazli/ScRNA_TIMP1_Bladder)

<img width="3600" height="1500" alt="QC63" src="https://github.com/user-attachments/assets/a57c3a57-d875-4de6-80a1-4755a1a39b2c" />
<img width="3600" height="2400" alt="QC186" src="https://github.com/user-attachments/assets/682a5bdf-5677-4f30-8831-61079de457c7" />
<img width="4200" height="3000" alt="QC105" src="https://github.com/user-attachments/assets/b6e935d2-18a9-4920-9a5c-6ee135cbee33" />

## Overview

This repository accompanies a study of **TIMP1 transcriptional heterogeneity in bladder carcinoma (BLCA)** using publicly available single-cell RNA sequencing (scRNA-seq) data. The manuscript also reports an **independent Egyptian female case–control genetic association study** of TIMP1 variants. These are separate cohorts and **no participant-level genotype–transcriptome linkage is available**.

The single-cell analyses characterize cell populations and TIMP1 expression patterns; they do not independently establish clinical biomarker performance, malignant identity for individual epithelial cells, functional stemness, or directional lineage transitions.

## Public input datasets

| Specimens | Source | Identifiers |
| --- | --- | --- |
| Tumor P1 | 10x Genomics, *40k Mixture of Dissociated Tumor Cells from 4 Donors, Multiplexed Samples, 4 Probe Barcodes* | [Official dataset page](https://www.10xgenomics.com/datasets/40k-mixture-of-dissociated-tumor-cells-from-4-donors-multiplexed-samples-4-probe-barcodes-1-standard); exact downloaded file(s) and probe barcode remain to be verified |
| Tumors P2–P3 | NCBI GEO | [GSE130001](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130001) |
| Tumors P4–P10; urine U1–U8 | NCBI GEO | [GSE267718](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267718) |
| Adjacent nonmalignant tissues N1–N3 | NCBI GEO | [GSE129845](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129845) |

The integrated analysis comprises **21 specimens** (10 tumors, 3 adjacent nonmalignant tissues, and 8 urine specimens), with **130,296 cells before quality control and 65,039 retained afterward**. Because specimen category and source study are partly or fully confounded, cross-category comparisons require caution.

Download original data from the primary repositories and consult their respective reuse conditions. This repository should not be assumed to contain redistributed raw sequencing files or human participant-level genotype records.

## P1 provenance and sample-level cautions

The reported source for tumor specimen **P1** is the 10x Genomics public dataset [*40k Mixture of Dissociated Tumor Cells from 4 Donors, Multiplexed Samples, 4 Probe Barcodes*](https://www.10xgenomics.com/datasets/40k-mixture-of-dissociated-tumor-cells-from-4-donors-multiplexed-samples-4-probe-barcodes-1-standard). The provider describes a multiplexed mixture of four dissociated tumor specimens, including a bladder transitional cell carcinoma specimen, profiled using Chromium Fixed RNA Profiling (Flex). The reported **8,371 P1 cells before QC** match the previously documented bladder specimen cell count; this numerical match supports the proposed mapping but **does not by itself verify the exact input files or probe-barcode demultiplexing**.

**Required provenance check before a reproducibility claim:** record the exact downloaded matrix/file names, 10x dataset release or download date, the bladder-specific probe barcode, any filtering/demultiplexing command, and ideally checksums. P1 must not be treated as the entire four-tumor mixture. Donor identities and sample pairing must not be inferred from the labels P1–P10, N1–N3, and U1–U8.

The independent Egyptian female genetic association cohort is **not** the source of the scRNA-seq specimens; do not imply participant-level matching between cohorts.

## Repository organization

The public repository currently exposes eleven module directories:

```text
ScRNA_TIMP1_Bladder/
├── Module_01/
├── Module_02/
├── Module_03/
├── Module_04/
├── Module_05/
├── Module_06/
├── Module_07/
├── Module_08/
├── Module_09/
├── Module_10/
├── Module_11/
└── README.md
```

**Module descriptions and exact execution order require verification against the scripts inside each directory.** Do not assume that numeric module order is the required run order. Each module should ultimately document its entry-point script, inputs, outputs, dependencies, and relationship to manuscript figures or tables.

## Analysis overview

The manuscript describes the following analytical components. Their availability and precise implementation in this repository should be checked against the committed code.

1. **Quality control:** filtering of cells and genes using original RNA count data.
2. **Normalization and integration:** Seurat v5 workflow including SCTransform and canonical-correlation-analysis-based integration.
3. **Clustering and annotation:** clustering followed by SingleR annotation using the Human Primary Cell Atlas reference, supplemented by lineage-associated markers.
4. **TIMP1 characterization:** cell-type and specimen-level transcript detection and expression summaries.
5. **Exploratory visualization:** MAGIC smoothing and PHATE embedding, interpreted as descriptive tools rather than independent validation of lineage transitions.
6. **Exploratory state ordering:** Monocle3, with inferred paths interpreted as computational state relationships rather than observed temporal progression.
7. **Annotation checks:** non-imputed marker panels, specimen-level cell counts, and SingleR score separation, where the corresponding scripts and outputs are supplied.

**Expression-data distinction:** raw UMI counts are appropriate for transcript detection and quality-control metrics; normalized but non-imputed expression supports standard transcriptional analyses; MAGIC-smoothed values must be labeled as imputed and should not be treated as directly observed counts. Verify the actual input assay for each script, particularly PHATE and Monocle3.

## Getting started

1. Install R and the dependencies declared in the relevant module scripts.
2. Obtain the public datasets from the links above.
3. Inspect each module's scripts for input file paths, object requirements, and package versions.
4. Run the verified workflow in its documented order, recording the R session information and output logs.
5. Compare regenerated cell counts, specimen identifiers, and figure outputs with the manuscript.

There is **not yet a verified one-command installation or execution procedure** documented here. Avoid presenting the repository as fully executable until module entry points and dependencies have been checked.

## Monocle3 pseudotime reproducibility

The audited `Module_11.R` invokes `order_cells(cds, reduction_method = "UMAP")` without explicitly passing `root_cells` or `root_pr_nodes`. **The root used for the reported pseudotime is currently undocumented.** Do not infer a root from the plot or claim that the existing direction of pseudotime is reproducible.

If the saved ordered `cell_data_set` object is available, inspect `cds@principal_graph_aux[["UMAP"]]$root_pr_nodes` and archive the recorded root together with the object and relevant code. If the original root cannot be recovered, select a biologically justified root, explicitly record the corresponding cell IDs or principal node(s), rerun `order_cells()`, and regenerate affected plots and text. A new root is a **new analysis**, not a reconstruction of the original choice. Computational pseudotime does not establish observed lineage transitions.

## Publication release and environment status

**Status: code shared; complete end-to-end reproduction not yet verified.** The repository should include `sessionInfo.txt`, `package_versions.csv`, and preferably `renv.lock` (plus pinned Python dependencies where relevant). Environment files generated today document the **current** environment, not necessarily the historical environment used for the manuscript.

After verifying the scripts, inputs, outputs, and figure mapping, create a tagged GitHub release (for example `v1.0.0`) and record the **actual tag, commit SHA, release URL, and date** here and in the manuscript. No release tag or archival DOI is asserted in this README until one is published.

| Reproducibility item | Current status |
| --- | --- |
| P1 dataset landing page | Documented; exact file and barcode pending verification |
| GEO study accessions | Documented; specimen-level GSM/file mapping pending |
| R package versions | Capture script supplied; verified historical versions unavailable |
| Monocle3 root | Not documented in inspected script; recovery or rerun required |
| Verified module execution order | Pending |
| Publication tag / release | Pending |

## Interpretation and limitations

Cell-type annotations are computational predictions. In particular, epithelial identity does not itself establish malignancy; uncommon stem-like, neuronal-like, and progenitor-like populations require additional validation. SingleR `delta.next` reflects separation of reference-label scores, **not a calibrated probability**. MAGIC, PHATE, and Monocle3 do not demonstrate temporal direction, lineage tracing, tumor progression, or treatment response. The genetic case–control findings are independent of the scRNA-seq cohort and do not demonstrate that the studied variants alter cell-specific TIMP1 expression.

## Data and code availability

Public scRNA-seq input data are available from GEO accessions **GSE130001, GSE267718, and GSE129845**, and the [10x Genomics four-donor multiplexed tumor dataset](https://www.10xgenomics.com/datasets/40k-mixture-of-dissociated-tumor-cells-from-4-donors-multiplexed-samples-4-probe-barcodes-1-standard) used as the reported source for P1 (exact input files and probe-barcode assignment still require verification). Study-associated code is hosted in this repository. The availability of particular scripts, processed data, and results should be verified from the current repository contents.

## Citation

Please cite the associated manuscript when published. Add the final title, authors, journal, year, DOI, and a versioned repository release or archival DOI here.

## Contact

For questions, use the repository's GitHub Issues page, if enabled.

## License

No license is asserted by this README. Add a `LICENSE` file if redistribution and reuse permissions are intended.

### Specimen-to-accession mapping

| Specimen | Category | Public source |
|---|---|---|
| P1 | Tumor | [10x Genomics multiplexed four-donor tumor dataset](https://www.10xgenomics.com/datasets/40k-mixture-of-dissociated-tumor-cells-from-4-donors-multiplexed-samples-4-probe-barcodes-1-standard); original input filename and probe barcode not yet verified |
| P2, P3 | Tumor | [GSE130001](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130001) |
| P4, P5, P6, P7, P8, P9, P10 | Tumor | [GSE267718](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267718) |
| N1, N2, N3 | Adjacent nonmalignant tissue | [GSE129845](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129845) |
| U1, U2, U3, U4, U5, U6, U7, U8 | Urine | [GSE267718](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267718) |

These are **study-level mappings**, not verified GEO sample (GSM) identifiers. A complete provenance table should additionally record GSM or file accession, donor identifier if known, library chemistry, pairing, and original filename. Do not infer donors from specimen labels.

### Analysis scripts and parameter audit

| Analysis | Script | Settings observed in code | Remaining reproducibility information |
|---|---|---|---|
| SCTransform and integration | `Module_05.R` | RNA layers split by `Batch`; `SCTransform(BLCA_M5)`; `RunPCA`; initial `RunUMAP(dims=1:30)`; `IntegrateLayers(method=CCAIntegration, normalization.method="SCT", verbose=FALSE)`; `FindNeighbors(reduction="integrated.dr",dims=1:30)`; `FindClusters(resolution=0.5)`; integrated UMAP/t-SNE using dimensions 1:30; memory limit 25 GiB | Exact Seurat/SeuratObject/sctransform versions, variable-feature selection, integration defaults, random seeds, batch definitions |
| SingleR | `Module_06.R` | `celldex::HumanPrimaryCellAtlasData()`; `SingleR(test=as.SingleCellExperiment(...),ref=ref,labels=ref$label.main)` and separate fine-label call with `ref$label.fine`; `plotScoreHeatmap(max.labels=25,clusters=seurat_clusters)` | Package/reference versions, exact assay used in converted object, final label reconciliation and confidence thresholds |
| MAGIC | `Module_07.R` | `magic(t(BLCA_M7@assays$RNA$data), knn=15, seed=42)`; imputed values subsequently stored in a Seurat assay for exploratory plots | Rmagic/Python MAGIC versions, all unspecified defaults, gene set passed to MAGIC, explicit record of downstream figures using smoothed values |
| PHATE | `Module_07.R` | `phate(t(BLCA_M7@assays$RNA$data))`; embedding extracted from `$embedding` | phateR/Python phate versions, default kNN/decay/t/seed and environment; confirm whether `BLCA_M7` at this point is the original non-imputed object (the code uses `BLCA_M7`, not `BLCA_MAGIC`) |
| Pathway module scoring | `Module_09.R` | `msigdbr(category="H")`; Hallmark `ANGIOGENESIS`, `APOPTOSIS`, and `INFLAMMATORY_RESPONSE` gene sets; `AddModuleScore(..., assay="SCT", features=list(gene_set), search=TRUE)` | msigdbr/MSigDB release, Seurat version, score control-feature settings/defaults and random seed; these are module scores, **not** GSEA enrichment p-values |
| Monocle3 | `Module_11.R` | `as.cell_data_set(BLCA_M11)`; `cluster_cells`; `learn_graph`; `order_cells(reduction_method="UMAP")`; pseudotime plotted and exported | Monocle3 version/commit, graph parameters, root cell/principal node, whether root was selected interactively, lineage-specific analysis and donor replication |

**Monocle3 caution:** `order_cells()` in the inspected script does not supply `root_cells` or `root_pr_nodes`. An interactive root choice may have occurred, but the root identity is **not recorded in the script**. Therefore, pseudotime direction is not reproducible from this script alone, and pooled paths must not be described as demonstrated lineage transitions.

**MAGIC/PHATE caution:** the PHATE call explicitly references `BLCA_M7@assays$RNA$data`, whereas the MAGIC output is assigned to a separate object `BLCA_MAGIC`. This suggests the shown PHATE call does not use MAGIC-smoothed input, but this should be checked against the saved object's state and actual execution. The use of direct `@assays$RNA$data` access may require adaptation for Seurat v5 layered assays.

**SingleR caution:** the code includes multiple annotation attempts and object names; it should be cleaned into one verified, noninteractive pipeline before claiming exact reproduction of the final cell labels.

### Capture exact software versions

Run the following in the same R environment used to regenerate the figures and commit the resulting text file:

```r
packages <- c("Seurat", "SeuratObject", "sctransform", "glmGamPoi",
              "SingleR", "celldex", "SingleCellExperiment", "Rmagic",
              "phateR", "reticulate", "msigdbr", "monocle3",
              "SeuratWrappers")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
versions <- vapply(packages, function(p) {
  if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p))
  else "NOT INSTALLED"
}, character(1))
write.csv(data.frame(package=names(versions), version=unname(versions)),
          "package_versions.csv", row.names=FALSE)
if (requireNamespace("reticulate", quietly=TRUE)) {
  writeLines(capture.output(reticulate::py_config()), "python_environment.txt")
}
```

"package","version"
"Seurat","5.4.0"
"SeuratObject","5.3.0"
"sctransform","0.4.3"
"glmGamPoi","1.18.0"
"SingleR","2.8.0"
"celldex","1.16.0"
"SingleCellExperiment","1.28.1"
"Rmagic","NOT INSTALLED"
"phateR","1.0.7"
"reticulate","1.45.0"
"msigdbr","25.1.1"
"monocle3","1.4.26"
"SeuratWrappers","0.4.0"

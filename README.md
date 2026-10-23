# ScRNA_TIMP1_Bladder

<img width="3600" height="1500" alt="QC63" src="https://github.com/user-attachments/assets/a57c3a57-d875-4de6-80a1-4755a1a39b2c" />
<img width="3600" height="2400" alt="QC186" src="https://github.com/user-attachments/assets/682a5bdf-5677-4f30-8831-61079de457c7" />
<img width="4200" height="3000" alt="QC105" src="https://github.com/user-attachments/assets/b6e935d2-18a9-4920-9a5c-6ee135cbee33" />

TIMP1 expression heterogeneity in bladder carcinoma
Repository: ScRNA_TIMP1_Bladder
Overview
This repository accompanies a study of TIMP1 transcriptional heterogeneity in bladder carcinoma (BLCA) using publicly available single-cell RNA sequencing (scRNA-seq) data. The manuscript also reports an independent Egyptian female case–control genetic association study of TIMP1 variants. These are separate cohorts and no participant-level genotype–transcriptome linkage is available.

The single-cell analyses characterize cell populations and TIMP1 expression patterns; they do not independently establish clinical biomarker performance, malignant identity for individual epithelial cells, functional stemness, or directional lineage transitions.
Public input datasets
Specimens	Source	Identifiers
Tumor P1	10x Genomics public dataset	Exact dataset URL/version: to be documented
Tumors P2–P3	NCBI GEO	GSE130001
Tumors P4–P10; urine U1–U8	NCBI GEO	GSE267718
Adjacent nonmalignant tissues N1–N3	NCBI GEO	GSE129845


The integrated analysis comprises 21 specimens (10 tumors, 3 adjacent nonmalignant tissues, and 8 urine specimens), with 130,296 cells before quality control and 65,039 retained afterward. Because specimen category and source study are partly or fully confounded, cross-category comparisons require caution.

Download original data from the primary repositories and consult their respective reuse conditions. This repository should not be assumed to contain redistributed raw sequencing files or human participant-level genotype records.
Repository organization
The public repository currently exposes eleven module directories:
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
Module descriptions and exact execution order require verification against the scripts inside each directory. Do not assume that numeric module order is the required run order. Each module should ultimately document its entry-point script, inputs, outputs, dependencies, and relationship to manuscript figures or tables.
Analysis overview
The manuscript describes the following analytical components. Their availability and precise implementation in this repository should be checked against the committed code.
1. Quality control: filtering of cells and genes using original RNA count data.
2. Normalization and integration: Seurat v5 workflow including SCTransform and canonical-correlation-analysis-based integration.
3. Clustering and annotation: clustering followed by SingleR annotation using the Human Primary Cell Atlas reference, supplemented by lineage-associated markers.
4. TIMP1 characterization: cell-type and specimen-level transcript detection and expression summaries.
5. Exploratory visualization: MAGIC smoothing and PHATE embedding, interpreted as descriptive tools rather than independent validation of lineage transitions.
6. Exploratory state ordering: Monocle3, with inferred paths interpreted as computational state relationships rather than observed temporal progression.
7. Annotation checks: non-imputed marker panels, specimen-level cell counts, and SingleR score separation, where the corresponding scripts and outputs are supplied.
Expression-data distinction: raw UMI counts are appropriate for transcript detection and quality-control metrics; normalized but non-imputed expression supports standard transcriptional analyses; MAGIC-smoothed values must be labeled as imputed and should not be treated as directly observed counts. Verify the actual input assay for each script, particularly PHATE and Monocle3.
Reproducibility checklist
Before considering the workflow fully reproducible, ensure the repository includes or documents:
- [ ] A source URL and version for the P1 10x Genomics dataset.
- [ ] A verified description and execution command for every module.
- [ ] R version, package versions, and installation instructions (sessionInfo() or renv.lock).
- [ ] Input paths, file formats, and required Seurat object metadata columns.
- [ ] Explicit quality-control thresholds, random seeds, and computational parameters.
- [ ] A map from scripts to manuscript and supplementary figures/tables.
- [ ] A distinction between analyses using raw, normalized non-imputed, and MAGIC-imputed expression.
- [ ] Instructions for regenerating intermediate objects and final outputs.
- [ ] A release tag or commit hash identifying the version used for publication.
Getting started
1. Install R and the dependencies declared in the relevant module scripts.
2. Obtain the public datasets from the links above.
3. Inspect each module's scripts for input file paths, object requirements, and package versions.
4. Run the verified workflow in its documented order, recording the R session information and output logs.
5. Compare regenerated cell counts, specimen identifiers, and figure outputs with the manuscript.
There is not yet a verified one-command installation or execution procedure documented here. Avoid presenting the repository as fully executable until module entry points and dependencies have been checked.
Interpretation and limitations
Cell-type annotations are computational predictions. In particular, epithelial identity does not itself establish malignancy; uncommon stem-like, neuronal-like, and progenitor-like populations require additional validation. SingleR delta.next reflects separation of reference-label scores, not a calibrated probability. MAGIC, PHATE, and Monocle3 do not demonstrate temporal direction, lineage tracing, tumor progression, or treatment response. The genetic case–control findings are independent of the scRNA-seq cohort and do not demonstrate that the studied variants alter cell-specific TIMP1 expression.
Data and code availability
Public scRNA-seq input data are available from GEO accessions GSE130001, GSE267718, and GSE129845, and a public 10x Genomics dataset whose exact source should be specified. Study-associated code is hosted in this repository. The availability of particular scripts, processed data, and results should be verified from the current repository contents.
Citation
Please cite the associated manuscript when published. Add the final title, authors, journal, year, DOI, and a versioned repository release or archival DOI here.
Contact
For questions, use the repository's GitHub Issues page, if enabled.
License
No license is asserted by this README. Add a LICENSE file if redistribution and reuse permissions are intended.

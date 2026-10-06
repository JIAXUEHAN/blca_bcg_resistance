# ==============================================================================
# 00_Reconstruct_scRNA_master.R
# BCG-NMIBC / Cancer Communications
#
# Purpose
# -------
# Reconstruct and document the analysis path that generated the current
# `seurat_filtered.rds` master object from the historical scripts:
#   - scRNA.R
#   - washdata.R
#   - Fig1.R
#
# IMPORTANT
# ---------
# 1) This script is a cleaned reconstruction of the historical workflow.
# 2) It starts from `merged.rds`, because the original scripts contain
#    placeholder/truncated filenames for the two upstream source RDS files.
# 3) `merged.rds` already represents the merge of:
#       - public NMIBC scRNA-seq dataset GSE269877
#       - four internal NMIBC scRNA-seq samples
# 4) The hard-coded cluster-to-cell-type maps below reproduce the historical
#    annotation logic. If clustering is re-run under different package versions,
#    marker expression MUST be re-checked before reusing cluster IDs verbatim.
# 5) The historical final QC branch used:
#       nFeature_RNA > 500
#       nFeature_RNA < 6000
#       nCount_RNA   > 1000
#    without an additional percent.mt cutoff.
#    Earlier exploratory code tested percent.mt < 15, but that was not the
#    contiguous branch used immediately before Harmony integration.
#
# Recommended use for the paper:
#   - Keep this script as the provenance/reconstruction script.
#   - Use the existing validated `seurat_filtered.rds` as the frozen master
#     object for downstream analyses.
# ==============================================================================


# ==============================================================================
# 0. Packages and paths
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(harmony)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

set.seed(1234)

# Run from:
# H:/韩佳雪/cancer communications/
project_dir <- "."

# Historical intermediate object:
merged_rds <- "../NMIBC_rds_standard/merged.rds"

# Output directories
dir.create("01_scRNA/00_objects", recursive = TRUE, showWarnings = FALSE)
dir.create("01_scRNA/01_QC", recursive = TRUE, showWarnings = FALSE)
dir.create("01_scRNA/02_annotation", recursive = TRUE, showWarnings = FALSE)
dir.create("09_Tables", recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 1. Load the historical merged object
# ==============================================================================
#
# Historical provenance:
#   public object: GSE269877
#   internal object: four internal samples
#
# The original scRNA.R removed reductions before merging and joined clinical
# metadata by `orig.ident`. Public cells were labeled dataset = "GSE269877";
# internal cells were subsequently labeled dataset = "Internal".
#
# In the merged object used by washdata.R:
#   - public counts were stored in assay `RNACleaned`
#   - internal four-sample counts were stored in assay `RNA`
# ==============================================================================

raw_obj <- readRDS(merged_rds)

stopifnot("RNA" %in% Assays(raw_obj))
stopifnot("RNACleaned" %in% Assays(raw_obj))

cat("Loaded merged object:\n")
print(dim(raw_obj))
print(table(raw_obj$dataset, useNA = "ifany"))


# ==============================================================================
# 2. Recover public and internal raw-count matrices
# ==============================================================================
#
# This follows the historical washdata.R logic:
#   Internal four samples -> RNA assay
#   Public GSE269877      -> RNACleaned assay
#
# Public and internal matrices are restricted to their common genes before
# creating a single clean RNA assay.
# ==============================================================================

extract_counts_from_assay5 <- function(assay_obj) {

  if (inherits(assay_obj, "Assay5")) {

    layer_names <- grep("^counts", names(assay_obj@layers), value = TRUE)

    if (length(layer_names) == 0) {
      stop("No counts layer found.")
    }

    feature_names <- rownames(assay_obj)

    mats <- lapply(layer_names, function(lyr) {

      mat <- assay_obj@layers[[lyr]]

      # Historical scripts restored row names from the assay feature list.
      if (is.null(rownames(mat)) || length(rownames(mat)) != nrow(mat)) {
        if (nrow(mat) > length(feature_names)) {
          stop("Layer contains more rows than the assay feature list.")
        }
        rownames(mat) <- feature_names[seq_len(nrow(mat))]
      }

      mat
    })

    common <- Reduce(intersect, lapply(mats, rownames))

    mats <- lapply(
      mats,
      function(m) m[common, , drop = FALSE]
    )

    out <- do.call(cbind, mats)

  } else {

    out <- GetAssayData(assay_obj, layer = "counts")
  }

  as(out, "dgCMatrix")
}


# Internal samples
local_mat <- extract_counts_from_assay5(raw_obj[["RNA"]])

# Public GSE269877
public_mat <- extract_counts_from_assay5(raw_obj[["RNACleaned"]])

cat("Public cells:", ncol(public_mat), "\n")
cat("Internal cells:", ncol(local_mat), "\n")

# Historical script reported approximately:
#   250,229 public cells + 25,857 internal cells = 276,086 total cells

final_common_genes <- intersect(
  rownames(public_mat),
  rownames(local_mat)
)

cat("Common genes:", length(final_common_genes), "\n")

combined_counts <- cbind(
  public_mat[final_common_genes, , drop = FALSE],
  local_mat[final_common_genes, , drop = FALSE]
)

combined_counts <- as(combined_counts, "dgCMatrix")


# ==============================================================================
# 3. Restore metadata and rebuild a clean Seurat object
# ==============================================================================
#
# The historical script forced matrix and metadata to a shared sequential cell
# naming scheme after confirming equal dimensions. Here we first try the safer
# barcode-based alignment; if that is impossible, we reproduce the historical
# sequential alignment only when dimensions match exactly.
# ==============================================================================

meta_data <- raw_obj@meta.data

if (
  !is.null(colnames(combined_counts)) &&
  all(colnames(combined_counts) %in% rownames(meta_data))
) {

  meta_data <- meta_data[colnames(combined_counts), , drop = FALSE]

} else {

  if (ncol(combined_counts) != nrow(meta_data)) {
    stop(
      "Counts and metadata dimensions differ: ",
      ncol(combined_counts), " cells in counts vs ",
      nrow(meta_data), " rows in metadata."
    )
  }

  warning(
    "Original barcodes could not be matched directly. ",
    "Reproducing historical sequential cell-ID alignment."
  )

  cell_ids <- paste0("Cell_", seq_len(ncol(combined_counts)))
  colnames(combined_counts) <- cell_ids
  rownames(meta_data) <- cell_ids
}

meta_data$dataset <- ifelse(
  is.na(meta_data$dataset),
  "Internal",
  meta_data$dataset
)

seu <- CreateSeuratObject(
  counts = combined_counts,
  meta.data = meta_data,
  project = "NMIBC_GSE269877_Internal"
)

rm(public_mat, local_mat, combined_counts)
gc()


# ==============================================================================
# 4. Convert Ensembl IDs to gene symbols when needed
# ==============================================================================
#
# Historical scripts stripped Ensembl version suffixes, mapped ENSEMBL -> SYMBOL
# using org.Hs.eg.db, retained unmapped Ensembl IDs, and made duplicated symbols
# unique.
# ==============================================================================

features_raw <- rownames(seu)

if (mean(grepl("^ENSG", features_raw)) > 0.5) {

  ensembl_clean <- sub("\\..*$", "", features_raw)

  gene_map <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = ensembl_clean,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )

  gene_symbols <- ifelse(
    is.na(gene_map),
    ensembl_clean,
    as.character(gene_map)
  )

  gene_symbols <- make.unique(gene_symbols)

  # Rebuild the object rather than directly modifying Assay5 feature names.
  counts_tmp <- LayerData(seu, assay = "RNA", layer = "counts")
  rownames(counts_tmp) <- gene_symbols

  meta_tmp <- seu@meta.data

  seu <- CreateSeuratObject(
    counts = counts_tmp,
    meta.data = meta_tmp,
    project = "NMIBC_GSE269877_Internal"
  )

  rm(counts_tmp, meta_tmp)
  gc()
}


# ==============================================================================
# 5. Cell-level QC
# ==============================================================================

seu[["percent.mt"]] <- PercentageFeatureSet(
  seu,
  pattern = "^MT-|^mt-"
)

seu[["percent.ribo"]] <- PercentageFeatureSet(
  seu,
  pattern = "^RP[SL]|^rp[sl]"
)

seu[["percent.hb"]] <- PercentageFeatureSet(
  seu,
  pattern = "^HB[A-Z]|^hb[a-z]"
)

qc_before <- data.frame(
  n_cells = ncol(seu),
  n_features = nrow(seu)
)

write.csv(
  qc_before,
  "01_scRNA/01_QC/qc_before_filtering.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------
# Historical FINAL QC branch
# ------------------------------------------------------------------
# Important:
# Earlier exploratory code also tested percent.mt < 15.
# The contiguous branch immediately preceding the final Harmony
# workflow used only the thresholds below.
# ------------------------------------------------------------------

seu_filtered <- subset(
  seu,
  subset =
    nFeature_RNA > 500 &
    nFeature_RNA < 6000 &
    nCount_RNA > 1000
)

qc_after <- data.frame(
  n_cells = ncol(seu_filtered),
  n_features = nrow(seu_filtered)
)

write.csv(
  qc_after,
  "01_scRNA/01_QC/qc_after_filtering.csv",
  row.names = FALSE
)

saveRDS(
  seu_filtered,
  "01_scRNA/00_objects/01_merged_cleaned.rds",
  compress = FALSE
)

rm(seu)
gc()


# ==============================================================================
# 6. Normalization, PCA, and Harmony integration
# ==============================================================================
#
# Historical parameters:
#   LogNormalize scale.factor = 10,000
#   2,000 variable genes (vst)
#   PCA = 30 PCs
#
# Several Harmony settings were explored.
# The last sample-level branch before the coarse annotation used orig.ident,
# theta = 4, lambda = 1.5, max_iter = 20.
# ==============================================================================

seu_filtered <- NormalizeData(
  seu_filtered,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

seu_filtered <- FindVariableFeatures(
  seu_filtered,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

seu_filtered <- ScaleData(
  seu_filtered,
  features = VariableFeatures(seu_filtered),
  verbose = FALSE
)

seu_filtered <- RunPCA(
  seu_filtered,
  npcs = 30,
  verbose = FALSE
)

set.seed(1234)
seu_filtered <- RunHarmony(
  object = seu_filtered,
  group.by.vars = "orig.ident",
  theta = 4,
  lambda = 1.5,
  max_iter = 20,
  reduction.save = "harmony",
  verbose = FALSE
)

seu_filtered <- FindNeighbors(
  seu_filtered,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

# Historical coarse-clustering solution
seu_filtered <- FindClusters(
  seu_filtered,
  resolution = 0.2,
  verbose = FALSE
)

set.seed(1234)
seu_filtered <- RunUMAP(
  seu_filtered,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap",
  verbose = FALSE
)


# ==============================================================================
# 7. First-pass major cell-type annotation
# ==============================================================================
#
# Historical 23-cluster map from washdata.R.
# These labels were used to construct `merged_harmony_annotated.rds`, which was
# subsequently passed to scDblFinder.
# ==============================================================================

main_type_map <- c(
  "0"  = "Epithelial",
  "1"  = "Epithelial",
  "2"  = "Epithelial",
  "3"  = "T/NK cells",
  "4"  = "Epithelial",
  "5"  = "Epithelial",
  "6"  = "Myeloid/Macrophage",
  "7"  = "Fibroblasts",
  "8"  = "Epithelial",
  "9"  = "Epithelial",
  "10" = "Endothelial",
  "11" = "Epithelial",
  "12" = "Epithelial",
  "13" = "Plasma cells",
  "14" = "Pericytes/SMC",
  "15" = "Epithelial",
  "16" = "Epithelial",
  "17" = "Epithelial",
  "18" = "Plasma cells",
  "19" = "Epithelial",
  "20" = "Epithelial",
  "21" = "Epithelial",
  "22" = "Epithelial"
)

cluster_ids_now <- sort(unique(as.character(seu_filtered$seurat_clusters)))

if (!all(cluster_ids_now %in% names(main_type_map))) {
  warning(
    "Current cluster IDs do not exactly match the historical coarse map. ",
    "Inspect markers before applying annotation."
  )
}

seu_filtered$cell_type_main <- unname(
  main_type_map[as.character(seu_filtered$seurat_clusters)]
)

saveRDS(
  seu_filtered,
  "01_scRNA/00_objects/02_merged_harmony_annotated.rds",
  compress = FALSE
)


# ==============================================================================
# 8. Doublet removal with scDblFinder
# ==============================================================================
#
# Historical workflow:
#   convert Seurat -> SingleCellExperiment
#   scDblFinder(..., samples = "orig.ident")
#   retain only scDblFinder.class == "singlet"
# ==============================================================================

sce <- as.SingleCellExperiment(seu_filtered)

set.seed(1234)
sce <- scDblFinder(
  sce,
  samples = "orig.ident"
)

seu_filtered$scDblFinder.class <- sce$scDblFinder.class

merged_harmony_clean <- subset(
  seu_filtered,
  subset = scDblFinder.class == "singlet"
)

cat(
  "Cells remaining after scDblFinder:",
  ncol(merged_harmony_clean),
  "\n"
)

rm(sce, seu_filtered)
gc()


# ==============================================================================
# 9. Re-cluster singlets at resolution 0.8
# ==============================================================================
#
# Historical Fig1.R:
#   seurat_raw <- merged_harmony_clean
#   neighbors: harmony PCs 1:30
#   clustering resolution: 0.8
#   UMAP: harmony PCs 1:30
# ==============================================================================

seurat_raw <- merged_harmony_clean

seurat_raw <- FindNeighbors(
  seurat_raw,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

seurat_raw <- FindClusters(
  seurat_raw,
  resolution = 0.8,
  verbose = FALSE
)

set.seed(1234)
seurat_raw <- RunUMAP(
  seurat_raw,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "umap",
  verbose = FALSE
)


# ==============================================================================
# 10. Historical 48-cluster -> broad lineage map
# ==============================================================================
#
# IMPORTANT:
# These IDs are historical cluster IDs. They should not be blindly reused if
# the graph clustering changes under a different software environment.
# ==============================================================================

correct_cluster_ids <- c(
  "0" = "Epithelial", "1" = "Epithelial", "2" = "Epithelial", "3" = "Epithelial",
  "4" = "Epithelial", "5" = "Epithelial", "6" = "T/NK cells", "7" = "Epithelial",
  "8" = "Epithelial", "9" = "Epithelial", "10" = "T/NK cells", "11" = "Myeloid",
  "12" = "Epithelial", "13" = "Myeloid", "14" = "Epithelial", "15" = "B cells",
  "16" = "Epithelial", "17" = "Epithelial", "18" = "Epithelial", "19" = "Endothelial",
  "20" = "Fibroblasts", "21" = "Epithelial", "22" = "Epithelial", "23" = "T/NK cells",
  "24" = "Epithelial", "25" = "T/NK cells", "26" = "Smooth Muscle Cells",
  "27" = "Epithelial", "28" = "Epithelial", "29" = "Plasma cells",
  "30" = "T/NK cells", "31" = "Epithelial", "32" = "Epithelial",
  "33" = "Epithelial", "34" = "Epithelial", "35" = "Myeloid",
  "36" = "Epithelial", "37" = "Epithelial", "38" = "Epithelial",
  "39" = "Epithelial", "40" = "Epithelial", "41" = "Epithelial",
  "42" = "Epithelial", "43" = "Epithelial", "44" = "Epithelial",
  "45" = "Plasma cells", "46" = "Epithelial", "47" = "Epithelial"
)

seurat_raw$cell_type_main_corrected <- unname(
  correct_cluster_ids[as.character(seurat_raw$seurat_clusters)]
)


# ==============================================================================
# 11. Sub-cluster each broad lineage
# ==============================================================================
#
# Historical subcluster workflow:
#   NormalizeData
#   2,000 HVGs
#   ScaleData
#   PCA 30
#   Harmony by Patient
#   neighbors / UMAP using harmony dims 1:20
#   resolution 0.3
#   FindAllMarkers(min.pct = 0.25, logfc.threshold = 0.25)
#
# The resulting cluster IDs were manually annotated using marker expression.
# ==============================================================================

detail_maps <- list(

  "B cells" = c(
    "0"  = "Mature B cells",
    "1"  = "Epithelial_contam",
    "2"  = "Mast_contam",
    "3"  = "T_NK_contam",
    "4"  = "LAMP3_DC",
    "5"  = "SPP1_Mac",
    "6"  = "Langerhans_cDC2",
    "7"  = "pDC",
    "8"  = "Plasma cells",
    "9"  = "FCN1_Mono",
    "10" = "cDC1"
  ),

  "Endothelial" = c(
    "0"  = "ACKR1_Vein_EC",
    "1"  = "SELE_Activated_EC",
    "2"  = "Tip_Angiogenic_EC",
    "3"  = "Arterial_EC",
    "4"  = "Epithelial_contam",
    "5"  = "Epithelial_contam",
    "6"  = "LowQuality_contam",
    "7"  = "Immune_contam",
    "8"  = "Lymphatic_EC",
    "9"  = "Proliferating_EC",
    "10" = "T_NK_contam"
  ),

  "Epithelial" = c(
    "0"  = "Inflammatory_Epi",
    "1"  = "Hypoxic_Epi",
    "2"  = "Differentiated_Epi",
    "3"  = "Basal_EMT_Epi",
    "4"  = "Secretory_Epi",
    "5"  = "Proliferating_Epi",
    "6"  = "Metabolic_Epi",
    "7"  = "S100A8_Inflam_Epi",
    "8"  = "Plasma_contam",
    "9"  = "LowQuality_contam",
    "10" = "Proliferating_Epi_2",
    "11" = "Cycling_Basal_Epi",
    "12" = "Specialized_Epi",
    "13" = "Early_Response_Epi"
  ),

  "Fibroblasts" = c(
    "0" = "Cartilage_like_Fib",
    "1" = "COL10A1_MyoCAF",
    "2" = "Matrix_Remodeling_Fib",
    "3" = "CXCL14_Secretory_Fib",
    "4" = "PI16_Progenitor_Fib",
    "5" = "Epithelial_contam"
  ),

  "Myeloid" = c(
    "0"  = "FOLR2_Macrophage",
    "1"  = "DUOX2_Inflammatory_Myeloid",
    "2"  = "SPP1_Macrophage",
    "3"  = "FCN1_Monocyte",
    "4"  = "cDC2",
    "5"  = "Cycling_Myeloid",
    "6"  = "Epithelial_contam",
    "7"  = "pDC",
    "8"  = "Low_quality_Doublets",
    "9"  = "Low_quality_lncRNA",
    "10" = "Fibroblast_contam"
  ),

  "Plasma cells" = c(
    "0" = "Plasma_cells",
    "1" = "Erythroid_contam",
    "2" = "Fibroblast_contam",
    "3" = "Mast_cells_contam",
    "4" = "Epithelial_contam"
  ),

  "Smooth Muscle Cells" = c(
    "0" = "Inflammatory_SMC",
    "1" = "Low_quality_SMC",
    "2" = "Matrix_synthetic_SMC",
    "3" = "Epithelial_contam_1",
    "4" = "Fibroblast_like_SMC",
    "5" = "Epithelial_contam_2",
    "6" = "Neural_like_SMC"
  ),

  "T/NK cells" = c(
    "0"  = "CD4_Naive_Memory_T",
    "1"  = "CD8_Effector_Memory_T",
    "2"  = "Treg",
    "3"  = "CXCL13_T_cells",
    "4"  = "Stromal_contam",
    "5"  = "NK_gamma_delta_T",
    "6"  = "CD16_NK_cells",
    "7"  = "Proliferating_T_NK",
    "8"  = "Tfh_Exhausted_T",
    "9"  = "B_cell_contam",
    "10" = "Low_quality_lncRNA"
  )
)

seurat_raw$cell_type_detailed <- as.character(
  seurat_raw$cell_type_main_corrected
)

subcluster_objects <- list()

for (ct in names(detail_maps)) {

  message("Subclustering: ", ct)

  sub_obj <- subset(
    seurat_raw,
    subset = cell_type_main_corrected == ct
  )

  sub_obj <- NormalizeData(
    sub_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  sub_obj <- FindVariableFeatures(
    sub_obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )

  sub_obj <- ScaleData(
    sub_obj,
    features = VariableFeatures(sub_obj),
    verbose = FALSE
  )

  sub_obj <- RunPCA(
    sub_obj,
    npcs = 30,
    verbose = FALSE
  )

  set.seed(1234)
  sub_obj <- RunHarmony(
    sub_obj,
    group.by.vars = "Patient",
    reduction = "pca",
    reduction.save = "harmony",
    verbose = FALSE
  )

  sub_obj <- FindNeighbors(
    sub_obj,
    reduction = "harmony",
    dims = 1:20,
    verbose = FALSE
  )

  sub_obj <- FindClusters(
    sub_obj,
    resolution = 0.3,
    verbose = FALSE
  )

  set.seed(1234)
  sub_obj <- RunUMAP(
    sub_obj,
    reduction = "harmony",
    dims = 1:20,
    verbose = FALSE
  )

  historical_map <- detail_maps[[ct]]
  cluster_now <- as.character(sub_obj$seurat_clusters)

  if (!all(unique(cluster_now) %in% names(historical_map))) {

    warning(
      "Subcluster IDs for ", ct,
      " differ from the historical annotation map. ",
      "Do NOT trust the hard-coded labels until markers are re-checked."
    )
  }

  sub_obj$historical_detail_label <- unname(
    historical_map[cluster_now]
  )

  # Write detailed labels back by barcode
  cells_to_update <- intersect(
    colnames(seurat_raw),
    colnames(sub_obj)
  )

  label_vec <- sub_obj$historical_detail_label
  names(label_vec) <- colnames(sub_obj)

  seurat_raw@meta.data[
    cells_to_update,
    "cell_type_detailed"
  ] <- label_vec[cells_to_update]

  subcluster_objects[[ct]] <- sub_obj
}

saveRDS(
  subcluster_objects,
  "01_scRNA/02_annotation/historical_subcluster_objects.rds",
  compress = FALSE
)


# ==============================================================================
# 12. Standardize detailed labels
# ==============================================================================
#
# This reproduces the historical detailed_map used in Fig1.R and explains the
# detailed labels present in the current master object.
# ==============================================================================

detailed_map <- c(
  "Mature B cells"            = "Mature B Cells",
  "Epithelial_contam"         = "SLITRK6+ Epithelial Cells",
  "Mast_contam"               = "TPSD1+ Mast Cells",
  "T_NK_contam"               = "GZMA+ T/NK Cells",
  "LAMP3_DC"                  = "LAMP3_DC",
  "SPP1_Mac"                  = "SPP1_Macrophage",
  "Langerhans_cDC2"           = "Langerhans_cDC2",
  "pDC"                       = "pDC",
  "Plasma cells"              = "Plasma Cells",
  "FCN1_Mono"                 = "FCN1_Monocyte",
  "cDC1"                      = "cDC1",

  "ACKR1_Vein_EC"             = "ACKR1_Vein_EC",
  "SELE_Activated_EC"         = "SELE_Activated_EC",
  "Tip_Angiogenic_EC"         = "Tip_Angiogenic_EC",
  "Arterial_EC"               = "Arterial_EC",
  "Immune_contam"             = "CCL21+ Lymphatic EC",
  "Lymphatic_EC"              = "Lymphatic_EC",
  "Proliferating_EC"          = "Proliferating_EC",

  "Inflammatory_Epi"          = "Inflammatory_Epi",
  "Hypoxic_Epi"               = "Hypoxic_Epi",
  "Differentiated_Epi"        = "Differentiated_Epi",
  "Basal_EMT_Epi"             = "Basal_EMT_Epi",
  "Secretory_Epi"             = "Secretory_Epi",
  "Proliferating_Epi"         = "Proliferating_Epi",
  "Metabolic_Epi"             = "Metabolic_Epi",
  "S100A8_Inflam_Epi"         = "S100A8_Inflam_Epi",
  "Plasma_contam"             = "IGHG1+ Plasma Cells",
  "Proliferating_Epi_2"       = "Proliferating_Epi_2",
  "Cycling_Basal_Epi"         = "Cycling_Basal_Epi",
  "Specialized_Epi"           = "Specialized_Epi",
  "Early_Response_Epi"        = "Early_Response_Epi",

  "Cartilage_like_Fib"        = "Cartilage_like_Fib",
  "COL10A1_MyoCAF"            = "COL10A1_MyoCAF",
  "Matrix_Remodeling_Fib"     = "Matrix_Remodeling_Fib",
  "CXCL14_Secretory_Fib"      = "CXCL14_Secretory_Fib",
  "PI16_Progenitor_Fib"       = "PI16_Progenitor_Fib",

  "FOLR2_Macrophage"          = "FOLR2_Macrophage",
  "DUOX2_Inflammatory_Myeloid"= "DUOX2_Inflammatory_Myeloid",
  "SPP1_Macrophage"           = "SPP1_Macrophage",
  "FCN1_Monocyte"             = "FCN1_Monocyte",
  "cDC2"                      = "cDC2",
  "Cycling_Myeloid"           = "Cycling_Myeloid",
  "Fibroblast_contam"         = "COL1A1+ Fibroblasts",

  "Plasma_cells"              = "Plasma Cells",
  "Erythroid_contam"          = "HBB+ Erythroid Cells",
  "Mast_cells_contam"         = "TPSB2+ Mast Cells",

  "Inflammatory_SMC"          = "Inflammatory_SMC",
  "Matrix_synthetic_SMC"      = "Matrix_synthetic_SMC",
  "Epithelial_contam_1"       = "KRT8+ Epithelial Cells",
  "Fibroblast_like_SMC"       = "THY1+ Fibroblasts",
  "Epithelial_contam_2"       = "CLDN3+ Epithelial Cells",
  "Neural_like_SMC"           = "NCAM1+ Neural-like Cells",

  "CD4_Naive_Memory_T"        = "CD4_Naive_Memory_T",
  "CD8_Effector_Memory_T"     = "CD8_Effector_Memory_T",
  "Treg"                      = "Treg",
  "CXCL13_T_cells"            = "CXCL13_T_cells",
  "Stromal_contam"            = "FAT1+ Stromal Cells",
  "NK_gamma_delta_T"          = "NK_gamma_delta_T",
  "CD16_NK_cells"             = "CD16_NK_cells",
  "Proliferating_T_NK"        = "Proliferating_T_NK",
  "Tfh_Exhausted_T"           = "Tfh_Exhausted_T",
  "B_cell_contam"             = "MS4A1+ B Cells",

  "LowQuality_contam"         = "LowQuality_Cells",
  "Low_quality_Doublets"      = "LowQuality_Doublets",
  "Low_quality_lncRNA"        = "LowQuality_lncRNA",
  "Low_quality_SMC"           = "LowQuality_SMC"
)

curr_detailed <- as.character(seurat_raw$cell_type_detailed)

seurat_raw$cell_type_detailed <- ifelse(
  curr_detailed %in% names(detailed_map),
  unname(detailed_map[curr_detailed]),
  curr_detailed
)


# ==============================================================================
# 13. Logical reassignment of broad cell types
# ==============================================================================

seurat_raw$cell_type_main_corrected <- case_when(

  grepl(
    "LowQuality|Low_quality",
    seurat_raw$cell_type_detailed,
    ignore.case = TRUE
  ) ~ "Low Quality Cells",

  grepl(
    "Mast",
    seurat_raw$cell_type_detailed,
    ignore.case = TRUE
  ) ~ "Mast Cells",

  grepl(
    "Plasma",
    seurat_raw$cell_type_detailed,
    ignore.case = TRUE
  ) ~ "Plasma Cells",

  grepl(
    "Erythroid",
    seurat_raw$cell_type_detailed,
    ignore.case = TRUE
  ) ~ "Erythroid Cells",

  grepl(
    "B Cells|Mature B",
    seurat_raw$cell_type_detailed,
    ignore.case = TRUE
  ) ~ "B Cells",

  grepl(
    "Macrophage|Monocyte|cDC|pDC|DC|Myeloid",
    seurat_raw$cell_type_detailed,
    ignore.case = TRUE
  ) ~ "Myeloid Cells",

  grepl(
    "T_cells|Treg|Memory_T|NK|T_NK|Tfh|Exhausted",
    seurat_raw$cell_type_detailed,
    ignore.case = TRUE
  ) ~ "T/NK Cells",

  grepl(
    "Epi|Epithelial",
    seurat_raw$cell_type_detailed,
    ignore.case = TRUE
  ) ~ "Epithelial Cells",

  grepl(
    "Fib|CAF|Stromal",
    seurat_raw$cell_type_detailed,
    ignore.case = TRUE
  ) ~ "Fibroblasts",

  grepl(
    "EC|Lymphatic",
    seurat_raw$cell_type_detailed,
    ignore.case = TRUE
  ) ~ "Endothelial Cells",

  grepl(
    "SMC",
    seurat_raw$cell_type_detailed,
    ignore.case = TRUE
  ) ~ "Smooth Muscle Cells",

  grepl(
    "Neural",
    seurat_raw$cell_type_detailed,
    ignore.case = TRUE
  ) ~ "Neural-like Cells",

  TRUE ~ "Other"
)


# ==============================================================================
# 14. Historical cleanup before final master object
# ==============================================================================

# Remove low-quality detailed groups
seurat_raw <- subset(
  seurat_raw,
  subset = cell_type_main_corrected != "Low Quality Cells"
)

# Historical figure workflow subsequently merged mast cells into myeloid
seurat_raw$cell_type_main_corrected <- ifelse(
  seurat_raw$cell_type_main_corrected == "Mast Cells",
  "Myeloid Cells",
  seurat_raw$cell_type_main_corrected
)

# Remove erythroid and neural-like groups before the final seurat_filtered object
keep_cells <- colnames(seurat_raw)[
  !seurat_raw$cell_type_main_corrected %in%
    c("Erythroid Cells", "Neural-like Cells")
]

seurat_filtered <- subset(
  seurat_raw,
  cells = keep_cells
)

Idents(seurat_filtered) <- "cell_type_main_corrected"


# ==============================================================================
# 15. Final reductions stored in the historical seurat_filtered.rds
# ==============================================================================
#
# The current master object contains:
#   pca
#   harmony
#   umap
#   tsne (historical object)
#   umap_uncorrected
#   harmony_soft
#   umap_harmony_soft
#   harmony_ident
#   umap_harmony_ident
#
# The following reproduces the three additional UMAP/Harmony reductions that
# were explicitly created immediately before `seurat_filtered.rds` was saved.
# ==============================================================================

# Uncorrected PCA UMAP
set.seed(1234)
seurat_filtered <- RunUMAP(
  seurat_filtered,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap_uncorrected",
  verbose = FALSE
)

# Soft Harmony by dataset
set.seed(1234)
seurat_filtered <- RunHarmony(
  seurat_filtered,
  group.by.vars = "dataset",
  theta = 0.5,
  reduction.save = "harmony_soft",
  verbose = FALSE
)

set.seed(1234)
seurat_filtered <- RunUMAP(
  seurat_filtered,
  reduction = "harmony_soft",
  dims = 1:30,
  reduction.name = "umap_harmony_soft",
  verbose = FALSE
)

# Harmony by sample identity
set.seed(1234)
seurat_filtered <- RunHarmony(
  seurat_filtered,
  group.by.vars = "orig.ident",
  theta = 1.0,
  reduction.save = "harmony_ident",
  verbose = FALSE
)

set.seed(1234)
seurat_filtered <- RunUMAP(
  seurat_filtered,
  reduction = "harmony_ident",
  dims = 1:30,
  reduction.name = "umap_harmony_ident",
  verbose = FALSE
)


# ==============================================================================
# 16. Final checks
# ==============================================================================

cat("\nFinal master dimensions:\n")
print(dim(seurat_filtered))

cat("\nDataset composition:\n")
print(table(seurat_filtered$dataset, useNA = "ifany"))

cat("\nTumor groups:\n")
print(table(seurat_filtered$Tumor.Category, useNA = "ifany"))

cat("\nMajor cell types:\n")
print(table(seurat_filtered$cell_type_main_corrected, useNA = "ifany"))

cat("\nDetailed cell types:\n")
print(table(seurat_filtered$cell_type_detailed, useNA = "ifany"))

cat("\nReductions:\n")
print(Reductions(seurat_filtered))


# Expected current frozen master:
#   13,972 features
#   223,229 cells
#
# Major cell-type counts in the currently validated object:
#   B Cells               1,075
#   Endothelial Cells       960
#   Epithelial Cells    183,091
#   Fibroblasts            2,460
#   Myeloid Cells          9,250
#   Plasma Cells           7,254
#   Smooth Muscle Cells      264
#   T/NK Cells            18,875
#
# If a complete rerun does not reproduce these counts, do not overwrite the
# frozen master. Investigate software-version/random-seed/cluster-ID differences.


# ==============================================================================
# 17. Save reconstructed master
# ==============================================================================

saveRDS(
  seurat_filtered,
  "01_scRNA/00_objects/seurat_filtered_reconstructed.rds",
  compress = FALSE
)

# Save sample-level summary for the manuscript
sample_summary <- seurat_filtered@meta.data %>%
  group_by(
    orig.ident,
    Patient,
    Tumor.Category,
    Groups,
    dataset,
    Grade,
    Stage,
    Sex
  ) %>%
  summarise(
    n_cells = n(),
    .groups = "drop"
  )

write.csv(
  sample_summary,
  "09_Tables/scRNA_master_sample_summary.csv",
  row.names = FALSE
)

# Freeze package/session information
writeLines(
  capture.output(sessionInfo()),
  "01_scRNA/00_objects/sessionInfo_scRNA_master.txt"
)

cat("\nDONE.\n")
cat("Reconstructed master saved to:\n")
cat("01_scRNA/00_objects/seurat_filtered_reconstructed.rds\n")

# ------------------------------------------------------------
# scRNAseq_seurat_pipeline.R
# A parameterized Seurat pipeline
# ------------------------------------------------------------
# What it does:
#  - Loads a Seurat object (RDS) OR a 10x Genomics folder
#  - QC (mito%, nCount_RNA, nFeature_RNA)
#  - Normalization (SCTransform OR LogNormalize)
#  - PCA -> neighbors -> UMAP -> clustering
#  - Marker discovery per cluster
#  - Saves figures and tables to ./outputs
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

set.seed(42)
message("Starting scRNA-seq Seurat pipeline...")

# -----------------------------
# User parameters
# -----------------------------
INPUT_PATH        <- "data/raw/example.rds"   # path to an RDS file or a 10x folder
IS_RDS            <- TRUE                     # TRUE if INPUT_PATH is a .rds file, FALSE for 10x folder
SPECIES           <- "human"                 # "human" or "mouse" (for mito gene pattern)
NORM_METHOD       <- "sct"                   # "sct" or "lognorm"
UMAP_DIMS         <- 1:30
CLUSTER_RES       <- 0.8

# QC thresholds
NFEATURE_LOW      <- 200
NFEATURE_HIGH     <- 6000
MT_HIGH           <- 15

# Output directories
OUT_DIR           <- "outputs"
FIG_DIR           <- file.path(OUT_DIR, "figures")
TAB_DIR           <- file.path(OUT_DIR, "tables")
RDS_DIR           <- "data/processed"

# -----------------------------
# Parameter validation
# -----------------------------
SPECIES <- tolower(SPECIES)
NORM_METHOD <- tolower(NORM_METHOD)

if (!SPECIES %in% c("human", "mouse")) {
  stop("SPECIES must be 'human' or 'mouse'")
}

if (!NORM_METHOD %in% c("sct", "lognorm")) {
  stop("NORM_METHOD must be 'sct' or 'lognorm'")
}

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper: save a ggplot with error handling
# -----------------------------
save_plot <- function(p, filename, width = 7, height = 5, dpi = 300) {
  out <- file.path(FIG_DIR, paste0(filename, ".png"))
  tryCatch({
    ggsave(out, plot = p, width = width, height = height, dpi = dpi)
    message("Saved: ", out)
  }, error = function(e) {
    warning("Failed to save plot: ", out, " | ", e$message)
  })
}

# -----------------------------
# Load data
# -----------------------------
message("Loading data...")
if (IS_RDS) {
  if (!file.exists(INPUT_PATH)) {
    stop("RDS file not found: ", INPUT_PATH)
  }
  seu <- readRDS(INPUT_PATH)
} else {
  if (!dir.exists(INPUT_PATH)) {
    stop("10x folder not found: ", INPUT_PATH)
  }
  mtx <- Read10X(INPUT_PATH)
  seu <- CreateSeuratObject(mtx)
}

# -----------------------------
# Basic QC
# -----------------------------
message("Computing QC metrics...")
mito_pat <- if (SPECIES == "mouse") "^mt-" else "^MT-"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mito_pat)

seu <- subset(
  seu,
  subset = nFeature_RNA > NFEATURE_LOW &
           nFeature_RNA < NFEATURE_HIGH &
           percent.mt < MT_HIGH
)

if (ncol(seu) == 0) {
  stop("No cells remaining after QC filtering. Adjust QC thresholds.")
}

# QC plots
vln <- VlnPlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
save_plot(vln, "qc_violin")

feat1 <- FeaturePlot(
  seu,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)
save_plot(feat1, "qc_featureplot")

saveRDS(seu, file.path(RDS_DIR, "seurat_qc.rds"))
message("Saved: ", file.path(RDS_DIR, "seurat_qc.rds"))

# -----------------------------
# Normalization + DR + Clustering
# -----------------------------
message("Normalization + dimensionality reduction...")
if (NORM_METHOD == "sct") {
  seu <- SCTransform(seu, verbose = FALSE) |>
    RunPCA(verbose = FALSE) |>
    FindNeighbors(dims = UMAP_DIMS) |>
    FindClusters(resolution = CLUSTER_RES) |>
    RunUMAP(dims = UMAP_DIMS)
} else {
  seu <- NormalizeData(seu) |>
    FindVariableFeatures() |>
    ScaleData() |>
    RunPCA() |>
    FindNeighbors(dims = UMAP_DIMS) |>
    FindClusters(resolution = CLUSTER_RES) |>
    RunUMAP(dims = UMAP_DIMS)
}

# UMAP plots
p_umap_label <- DimPlot(seu, reduction = "umap", label = TRUE) +
  ggtitle("UMAP - clusters")

p_umap_group <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("UMAP - cluster groups")

save_plot(p_umap_label, "umap_clusters_label")
save_plot(p_umap_group, "umap_clusters_group")

saveRDS(seu, file.path(RDS_DIR, "seurat_clustered.rds"))
message("Saved: ", file.path(RDS_DIR, "seurat_clustered.rds"))

# -----------------------------
# Marker discovery
# -----------------------------
message("Finding markers (per cluster)...")
markers <- FindAllMarkers(
  seu,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
) |>
  arrange(cluster, desc(avg_log2FC))

write_csv(markers, file.path(TAB_DIR, "markers_per_cluster.csv"))
message("Saved: ", file.path(TAB_DIR, "markers_per_cluster.csv"))

# view + save. top marker per cluster
top_markers <- markers |>
  group_by(cluster) |>
  slice_max(order_by = avg_log2FC, n = 1)

write_csv(top_markers, file.path(TAB_DIR, "top_markers_per_cluster.csv"))
message("Saved: ", file.path(TAB_DIR, "top_markers_per_cluster.csv"))

print(top_markers |> select(cluster, gene, avg_log2FC), n = 50)

message("Done.")
# ------------------------------------------------------------
# scRNAseq_seurat_pipeline.R
# A clean, parameterized Seurat pipeline
# ------------------------------------------------------------
# What it does:
#  - Loads a Seurat object (RDS) OR a 10x Genomics folder
#  - QC (mito%, nCount_RNA, nFeature_RNA)
#  - Normalization (SCTransform OR LogNormalize)
#  - PCA → neighbors → UMAP → clustering
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

# -----------------------------
# User parameters (edit freely)
# -----------------------------
INPUT_PATH        <- "data/raw"        # path to an RDS file or a 10x folder
IS_RDS            <- TRUE              # TRUE if INPUT_PATH is a .rds file, FALSE for 10x folder
SPECIES           <- "human"           # "human" or "mouse" (for mito gene pattern)
NORM_METHOD       <- "sct"             # "sct" or "lognorm"
UMAP_DIMS         <- 1:30
CLUSTER_RES       <- 0.8
OUT_DIR           <- "outputs"
FIG_DIR           <- file.path(OUT_DIR, "figures")
TAB_DIR           <- file.path(OUT_DIR, "tables")
RDS_DIR           <- "data/processed"

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RDS_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper: save a ggplot safely
# -----------------------------
save_plot <- function(p, filename, width = 7, height = 5, dpi = 300) {
  out <- file.path(FIG_DIR, paste0(filename, ".png"))
  ggsave(out, plot = p, width = width, height = height, dpi = dpi)
  message("Saved: ", out)
}

# -----------------------------
# Load data
# -----------------------------
message("Loading data...")
if (IS_RDS) {
  stopifnot(file.exists(INPUT_PATH))
  seu <- readRDS(INPUT_PATH)
} else {
  stopifnot(dir.exists(INPUT_PATH))
  mtx <- Read10X(INPUT_PATH)
  seu <- CreateSeuratObject(mtx)
}

# -----------------------------
# Basic QC
# -----------------------------
message("Computing QC metrics...")
mito_pat <- if (tolower(SPECIES) == "mouse") "^mt-" else "^MT-"
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = mito_pat)

# Adjust thresholds to your dataset
nFeature_low  <- 200
nFeature_high <- 6000
mt_high       <- 15

seu <- subset(
  seu,
  subset = nFeature_RNA > nFeature_low &
           nFeature_RNA < nFeature_high &
           percent.mt < mt_high
)

# QC plots
vln <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
save_plot(vln, "qc_violin")

feat1 <- FeaturePlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
save_plot(feat1, "qc_featureplot")

saveRDS(seu, file.path(RDS_DIR, "seurat_qc.rds"))

# -----------------------------
# Normalization + DR + Clustering
# -----------------------------
message("Normalization + dimensionality reduction...")
if (tolower(NORM_METHOD) == "sct") {
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
    FindNeighbors(dims = 1:20) |>
    FindClusters(resolution = CLUSTER_RES) |>
    RunUMAP(dims = 1:20)
}

# UMAP plots
p_umap_label <- DimPlot(seu, reduction = "umap", label = TRUE) + ggtitle("UMAP — clusters")
p_umap_split <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters") + ggtitle("UMAP — cluster groups")
save_plot(p_umap_label, "umap_clusters_label")
save_plot(p_umap_split, "umap_clusters_group")

saveRDS(seu, file.path(RDS_DIR, "seurat_clustered.rds"))

# -----------------------------
# Marker discovery
# -----------------------------
message("Finding markers (per cluster)...")
markers <- FindAllMarkers(
  seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25
) |>
  arrange(cluster, desc(avg_log2FC))

write_csv(markers, file.path(TAB_DIR, "markers_per_cluster.csv"))

# Quick view: top marker per cluster (prints to console)
top_markers <- markers |> group_by(cluster) |> slice_max(order_by = avg_log2FC, n = 1)
print(top_markers |> select(cluster, gene, avg_log2FC), n = 50)

message("Done.")

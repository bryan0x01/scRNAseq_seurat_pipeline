# Example configuration for the Seurat pipeline

INPUT_PATH  <- "data/raw/example_dataset.rds"
IS_RDS      <- TRUE
SPECIES     <- "human"

# Normalization
NORM_METHOD <- "sct"

# Dimensionality reduction
UMAP_DIMS   <- 1:30

# Clustering resolution
CLUSTER_RES <- 0.8

# Output directory
OUT_DIR     <- "outputs"
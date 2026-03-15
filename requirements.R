# Install required packages for the scRNA-seq pipeline

required_packages <- c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "readr"
)

install_if_missing <- function(pkg){
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}

invisible(lapply(required_packages, install_if_missing))

cat("All required packages installed.\n")
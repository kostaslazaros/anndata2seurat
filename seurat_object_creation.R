#!/usr/bin/env Rscript
library(Seurat)
args <- commandArgs(trailingOnly = TRUE)

print("Hello from R")

print("Running script with parameters:")
print(paste("basic_path", args[1], sep = ":"))
print(paste("has_umap", args[2], sep = ":"))

basic_path <- args[1]
has_umap <- args[2]

metadata_file <- paste(basic_path, "metadata.csv", sep = "/")
coordinates_file <- paste(basic_path, "umap_coordinates.csv", sep = "/")
seurat_file <- paste(basic_path, "pdx_raw.rds", sep = "/")

print(paste("Reading raw data from path", basic_path, sep = ":"))
raw_data <- Read10X(basic_path)

print(paste("Reading metadata file", metadata_file, sep = ":"))
metadata <- read.csv(metadata_file, row.names = 1)

print("Creating seurat object")
seurat_object <- CreateSeuratObject(counts = raw_data, meta.data = metadata)

if (has_umap != "") {
  print(paste("Adding umap to seurat file from", coordinates_file, sep = ":"))
  umap_coords <- read.csv(coordinates_file, row.names = 1)
  seurat_object[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(umap_coords), key = "UMAP_")
  # DimPlot(seurat_object, group.by = "cell_type")
}

print(paste("Saving seurat as", seurat_file, sep = ":"))
saveRDS(seurat_object, file = seurat_file)

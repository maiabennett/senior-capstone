# Navigate to the relevant RDS file
filerds <- file.choose()

# Read RDS from that file
seurat <- readRDS(filerds)

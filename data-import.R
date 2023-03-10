# GSE 150050, four tissues
## Read in the counts data
counts.15 <- read.csv("C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/Data/GSE150050/STAR_raw_counts.csv", row.names=1)

## Read in the metadata
metadata.15 <- read.csv("C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/Data/GSE150050/GSE150050_metadata.csv", sep = ";", row.names=1)

## SQL script is used to narrow metadata.15 to colon samples
metadata15 <- as.data.frame(metadata.15)
metadata.15.colon <- sqldf("SELECT * FROM metadata15 WHERE TISSUE = 'COLON'")

## Metadata indicates all colon samples have distinct identifier 'GNI' in cell ID
counts15 <- as.data.frame(counts.15)
counts.15.colon <- counts15[,grepl("GNI", colnames(counts15))]

## Seurat object creation
seurat.15 <- CreateSeuratObject(counts = counts.15)

# GSE 185224, normal gut
## Read in the counts data
counts.18 <- Read10X_h5("C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/Data/GSE185224/GSE185224_Donor1_filtered_feature_bc_matrix.h5", use.names=TRUE)

## Seurat object creation from Gene Expression subset (not Antibody Capture subset)
seurat.18 <- CreateSeuratObject(counts = counts.18$`Gene Expression`)

# GSE 125527, IBD gut
# Intestinal immune cells
dir.12.int <- "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/Data/GSE125527/UMI/HealthyI/"
files.12.int <- list.files(path = dir.12.int, pattern = ".tsv.gz", full.names = TRUE)

### Read in and aggregate  samples
counts.12.int <- read.csv(files.12.int[1],sep="\t", row.names=1)
for(i in 2:length(files.12.int)){
  counts.12b.int <- read.csv(files.12.int[i],sep="\t", row.names=1)
  counts.12.int <- rbind(counts.12.int, counts.12b.int)
}
counts.12.int <- t(counts.12.int)

## Seurat objects creation
seurat.12.int <- CreateSeuratObject(counts = counts.12.int)


## PBMC counts data from directory and append
dir.12.healthy <- "C:/Users/maiabennett/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/Data/GSE125527/UMI/HealthyPBMC/"
dir.12.UC <- "C:/Users/maiabennett/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/Data/GSE125527/UMI/UCPBMC/"
files.12.pbmc.healthy <- list.files(path = dir.12.healthy, pattern = ".tsv.gz", full.names = TRUE)
files.12.pbmc.UC <- list.files(path = dir.12.UC, pattern = ".tsv.gz", full.names = TRUE)

### Read in and aggregate healthy control samples
counts.12.healthy <- read.csv(files.12.pbmc.healthy[1],sep="\t", row.names=1)
for(i in 2:length(files.12.pbmc.healthy)){
  counts.12b.healthy <- read.csv(files.12.pbmc.healthy[i],sep="\t", row.names=1)
  counts.12.healthy <- rbind(counts.12.healthy, counts.12b.healthy)
}
counts.12.pbmc.healthy <- t(counts.12.healthy)

### Make metadata for later analysis
metadata.12.healthy <- data.frame(x1= colnames(counts.12.pbmc.healthy), x2 = "healthy")
colnames(metadata.12.healthy) <- c("barcode", "condition")

### Read in and aggregate UC samples
counts.12.UC <- read.csv(files.12.pbmc.UC[1],sep="\t", row.names=1)
for(i in 2:length(files.12.pbmc.UC)){
  counts.12b.UC <- read.csv(files.12.pbmc.UC[i],sep="\t", row.names=1)
  counts.12.UC <- rbind(counts.12.UC, counts.12b.UC)
}
counts.12.pbmc.UC <- t(counts.12.UC)

### Make metadata for later analysis
metadata.12.UC <- data.frame(x1= colnames(counts.12.pbmc.UC), x2 = "UC")
colnames(metadata.12.UC) <- c("barcode", "condition")
metadata.12.UC <- column_to_rownames(metadata.12.UC, loc = 1)

## Seurat objects creation
seurat.12.healthy <- CreateSeuratObject(counts = counts.12.pbmc.healthy, metadata = metadata.12.healthy)
seurat.12.UC <- CreateSeuratObject(counts = counts.12.pbmc.UC, metadata = metadata.12.UC)

## Create combined object
counts.12 <- rbind(t(counts.12.pbmc.healthy), t(counts.12.pbmc.UC))
counts.12 <- t(counts.12)
metadata.12 <- rbind(metadata.12.healthy, metadata.12.UC)
metadata.12 <- column_to_rownames(metadata.12, loc = 1)
seurat.12 <- CreateSeuratObject(counts = counts.12)
seurat.12 <- AddMetaData(seurat.12, metadata.12, col.name = "condition")

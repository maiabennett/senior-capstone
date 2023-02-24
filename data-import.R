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
## Read in the PBMC counts data from directory and append
dir.12 <- "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/Data/GSE125527/UMI/"
files.12 <- list.files(path = dir.12, pattern = ".tsv.gz", full.names = TRUE)
files.12.pbmc <- files.12[grep("pBMC",files.12)]
counts.12 <- read.csv(files.12.pbmc[1],sep="\t", row.names=1)
for(i in 2:length(files.12.pbmc)){
  counts.12b <- read.csv(files.12.pbmc[i],sep="\t", row.names=1)
  counts.12 <- rbind(counts.12, counts.12b)
}
counts.12.pbmcs <- t(counts.12)

## Read in the metadata
metadata.12 <- read.csv("C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/Data/GSE125527/GSE125527_cell_metadata.csv", row.names=1)

## SQL script used to narrow metadata.12 to PBMC samples only
metadata12 <- as.data.frame(metadata.12)
metadata.12.pbmcs <- sqldf("SELECT * FROM metadata12 WHERE tissue_assignment= 'PBMC'")

## Seurat object creation
seurat.12 <- CreateSeuratObject(counts = counts.12.pbmcs, meta.data=metadata.12.pbmcs)

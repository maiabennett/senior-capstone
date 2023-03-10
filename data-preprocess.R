
# GSE 150050, four tissues

## QC and valid cell selection
seurat.15[["percent.mt"]] <- PercentageFeatureSet(seurat.15, pattern = "^MT-")
VlnPlot(seurat.15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat.15 <- subset(seurat.15, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

## Normalization of data
seurat.15 <- NormalizeData(seurat.15)

### Feature selection
seurat.15 <- FindVariableFeatures(seurat.15, selection.method = "vst", nfeatures = 2000)

### Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.15), 10)

### Plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat.15)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat15-varfeatures.png", width= 900, height=480)
plot1
plot2
dev.off()

## Scaling of data
all.genes <- rownames(seurat.15)
seurat.15 <- ScaleData(seurat.15, features = all.genes)

## Perform linear dimension reduction
seurat.15 <- RunPCA(seurat.15, features = VariableFeatures(object = seurat.15))
print(seurat.15[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(seurat.15, reduction = "pca")
DimHeatmap(seurat.15, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat.15, dims = 1:15, cells = 500, balanced = TRUE)
seurat.15 <- JackStraw(seurat.15, num.replicate = 100)
seurat.15 <- ScoreJackStraw(seurat.15, dims = 1:20)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat15-jackstraw.png", width= 700, height=480)
JackStrawPlot(seurat.15, dims = 1:15)
dev.off()

## Save object for future recall
saveRDS(seurat.15, file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat15-preprocessed.rds")

# GSE 185224, normal gut

## QC and valid cell selection
seurat.18[["percent.mt"]] <- PercentageFeatureSet(seurat.18, pattern = "^MT-")
VlnPlot(seurat.18, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat.18 <- subset(seurat.18, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

## Normalization of data
seurat.18 <- NormalizeData(seurat.18)

### Feature selection
seurat.18 <- FindVariableFeatures(seurat.18, selection.method = "vst", nfeatures = 2000)

### Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.18), 10)

### Plot variable features with and without labels
plot3 <- VariableFeaturePlot(seurat.18)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat18-varfeatures.png", width= 900, height=480)
plot3 + plot4
dev.off()

## Scaling of data
all.genes <- rownames(seurat.18)
seurat.18 <- ScaleData(seurat.18, features = all.genes)

## Perform linear dimension reduction
seurat.18 <- RunPCA(seurat.18, features = VariableFeatures(object = seurat.18))
print(seurat.18[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(seurat.18, reduction = "pca")
DimHeatmap(seurat.18, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat.18, dims = 1:15, cells = 500, balanced = TRUE)
seurat.18 <- JackStraw(seurat.18, num.replicate = 100)
seurat.18 <- ScoreJackStraw(seurat.18, dims = 1:20)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat18-jackstraw.png", width= 700, height=480)
JackStrawPlot(seurat.18, dims = 1:15)
dev.off()

## Save object for future recall
saveRDS(seurat.18, file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat18-preprocessed.rds")
# GSE 125527, IBD gut
# Intestinal immune cells
## QC and valid cell selection
seurat.12.int[["percent.mt"]] <- PercentageFeatureSet(seurat.12.int, pattern = "^MT-")
VlnPlot(seurat.12.int, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat.15 <- subset(seurat.12.int, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

## Normalization of data
seurat.12.int <- NormalizeData(seurat.12.int)

### Feature selection
seurat.12.int <- FindVariableFeatures(seurat.12.int, selection.method = "vst", nfeatures = 2000)

### Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.12.int), 10)

### Plot variable features with and without labels
print(plot1 <- VariableFeaturePlot(seurat.12.int))
print(plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE))

## Scaling of data
all.genes <- rownames(seurat.12.int)
seurat.12.int <- ScaleData(seurat.12.int, features = all.genes)

## Perform linear dimension reduction
seurat.12.int <- RunPCA(seurat.12.int, features = VariableFeatures(object = seurat.12.int))
print(seurat.12.int[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(seurat.12.int, reduction = "pca")
DimHeatmap(seurat.12.int, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat.12.int, dims = 1:15, cells = 500, balanced = TRUE)
seurat.12.int <- JackStraw(seurat.12.int, num.replicate = 100)
seurat.12.int <- ScoreJackStraw(seurat.12.int, dims = 1:20)
JackStrawPlot(seurat.12.int, dims = 1:15)

## Save object for future recall
saveRDS(seurat.12.int, file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat12int-preprocessed.rds")

# GSE125527

# Intestinal immune cells
## QC and valid cell selection
seurat.12.int[["percent.mt"]] <- PercentageFeatureSet(seurat.12.int, pattern = "^MT-")
VlnPlot(seurat.12.int, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat.12.int <- subset(seurat.12.int, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)

## Normalization of data
seurat.12.int <- NormalizeData(seurat.12.int)

### Feature selection
seurat.12.int <- FindVariableFeatures(seurat.12.int, selection.method = "vst", nfeatures = 2000)

### Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.12.int), 10)

### Plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat.12.int)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12int-varfeatures.png", width= 900, height=480)
plot1 + plot2
dev.off()

## Scaling of data
all.genes <- rownames(seurat.12.int)
seurat.12.int <- ScaleData(seurat.12.int, features = all.genes)

## Perform linear dimension reduction
seurat.12.int <- RunPCA(seurat.12.int, features = VariableFeatures(object = seurat.12.int))
print(seurat.12.int[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(seurat.12.int, reduction = "pca")
DimHeatmap(seurat.12.int, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat.12.int, dims = 1:15, cells = 500, balanced = TRUE)
seurat.12.int <- JackStraw(seurat.12.int, num.replicate = 100)
seurat.12.int <- ScoreJackStraw(seurat.12.int, dims = 1:20)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12int-jackstraw.png", width= 700, height=480)
JackStrawPlot(seurat.12.int, dims = 1:15)
dev.off()

## Save object for future recall
saveRDS(seurat.12.int, file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat12int-preprocessed.rds")


# PBMCs by condition
seurat.12.healthy[["percent.mt"]] <- PercentageFeatureSet(seurat.12.healthy, pattern = "^MT-")
VlnPlot(seurat.12.healthy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat.12.healthy <- subset(seurat.12.healthy, subset = nFeature_RNA > 200 & nFeature_RNA < 4000)

seurat.12.UC[["percent.mt"]] <- PercentageFeatureSet(seurat.12.UC, pattern = "^MT-")
VlnPlot(seurat.12.UC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat.12.UC <- subset(seurat.12.UC, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

## Normalization of data
seurat.12.healthy <- NormalizeData(seurat.12.healthy)
seurat.12.UC <- NormalizeData(seurat.12.UC)

### Feature selection
seurat.12.healthy <- FindVariableFeatures(seurat.12.healthy, selection.method = "vst", nfeatures = 2000)
seurat.12.UC <- FindVariableFeatures(seurat.12.UC, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(seurat.12.healthy), 10)
top10 <- head(VariableFeatures(seurat.12.UC), 10)

### Plot individual variable features with and without labels
print(plot5 <- VariableFeaturePlot(seurat.12.healthy))
print(plot6 <- LabelPoints(plot = plot5, points = top10, repel = TRUE))
print(plot7 <- VariableFeaturePlot(seurat.12.UC))
print(plot8 <- LabelPoints(plot = plot7, points = top10, repel = TRUE))
png(file = "C:/Users/maiabennett/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/12healthy-varfeatures.png", width= 900, height=480)
plot5 + plot6
dev.off()
png(file = "C:/Users/maiabennett/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12UC-varfeatures.png", width= 900, height=480)
plot7 + plot8
dev.off()

### Normalize and identify variable features in combined Seurat object
seurat.12.list <- SplitObject(seurat.12, split.by = "condition")
seurat.12.list <- lapply(X = seurat.12.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures=2000)
})

### Select features that are repeatedly variable across datasets for integration
seurat.12.features <- SelectIntegrationFeatures(object.list = seurat.12.list)
seurat.12.anchors <- FindIntegrationAnchors(object.list = seurat.12.list, anchor.features = seurat.12.features)

### Create an 'integrated' data assay
seurat.12 <- IntegrateData(anchorset = seurat.12.anchors)
DefaultAssay(seurat.12) <- "integrated"

### Scaling and linear dimension reduction
seurat.12 <- ScaleData(seurat.12, verbose = FALSE)
seurat.12 <- RunPCA(seurat.12, npcs = 30, verbose = FALSE)
print(seurat.12[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(seurat.12, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat.12, dims = 1:15, cells = 500, balanced = TRUE)
seurat.12 <- JackStraw(seurat.12, num.replicate = 100)
seurat.12 <- ScoreJackStraw(seurat.12, dims = 1:20)
png(file = "C:/Users/maiabennett/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12-jackstraw.png", width= 700, height=480)
JackStrawPlot(seurat.12, dims = 1:15)
dev.off()

seurat.12.healthy <- ScaleData(seurat.12.healthy, verbose = FALSE)
seurat.12.healthy <- RunPCA(seurat.12.healthy, npcs = 30, verbose = FALSE)
print(seurat.12.healthy[["pca"]], dims = 1:5, nfeatures = 5)

seurat.12.UC <- ScaleData(seurat.12.UC, verbose = FALSE)
seurat.12.UC <- RunPCA(seurat.12.UC, npcs = 30, verbose = FALSE)
print(seurat.12.UC[["pca"]], dims = 1:5, nfeatures = 5)

## Save object for future recall
saveRDS(seurat.12.healthy, file = "C:/Users/maiabennett/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat12healthy-preprocessed.rds")
## Save object for future recall
saveRDS(seurat.12.UC, file = "C:/Users/maiabennett/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat12UC-preprocessed.rds")
## Save object for future recall
saveRDS(seurat.12, file = "C:/Users/maiabennett/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat12-preprocessed.rds")

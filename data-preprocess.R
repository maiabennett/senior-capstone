# GSE 150050, four tissues
## QC and valid cell selection
seurat.15[["percent.mt"]] <- PercentageFeatureSet(seurat.15, pattern = "^MT-")
VlnPlot(seurat.15, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seurat.15, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.15, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
seurat.15 <- subset(seurat.15, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## Normalization of data
seurat.15 <- NormalizeData(seurat.15)

### Feature selection
seurat.15 <- FindVariableFeatures(seurat.15, selection.method = "vst", nfeatures = 2000)

### Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat.15), 10)

### plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat.15)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## Scaling of data
all.genes <- rownames(seurat.15)
seurat.15 <- ScaleData(seurat.15, features = all.genes)

# GSE 185224, normal gut

# GSE 125527, IBD gut
# GSE 125527
# PBMCs by condition

## Cluster the cells
seurat.12 <- FindNeighbors(seurat.12, dims = 1:10)
seurat.12 <- FindClusters(seurat.12, resolution = 0.5)

### Look at cluster IDs of the first 5 cells
head(Idents(seurat.12), 5)

## Run non-linear dimensional reduction (UMAP/tSNE)
seurat.12 <- RunUMAP(seurat.12, dims = 1:10)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12-initialUMAP.png")
DimPlot(seurat.12, reduction = "umap", label = TRUE)
dev.off()
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12-initialUMAPcond.png")
DimPlot(seurat.12, reduction = "umap", group.by = "condition")
dev.off()

## Find markers that differentiate clusters from one another (to help assign identities)
DefaultAssay(seurat.12) <- "RNA"
seurat.12.markers <- FindAllMarkers(seurat.12, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.12 <- seurat.12.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.table(manual.curate.12, file="manual.curate.12.txt", row.names=TRUE, col.names=TRUE)

## Visualize specific markers to identify cell types
### ILC1s: CD127/IL7R+, CD161/KLRB1+, ID2+, T-bet/TBX21+, CD294/PTGDR2-, CD117/KIT-
## ILC2s: CD127/IL7R+, CD161/KLRB1+, ID2+, CD294/PTGDR2+, GATA3+, CD117/KIT+-
### ILC3s: CD127/IL7R+, CD161/KLRB1+, ID2+, CD294/PTGDR-, CD117/KIT+, IL23R+
### NKs: CD56/NCAM1+, EOMES+
### Ts: CD3/CD3D+ (definite), TBX21+
### Monocytes: CD14+
### Neutrophils, Basophils, Eosinophils, DCs: CD11b+
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12-initialmarkers.png", width= 700, height=480)
FeaturePlot(seurat.12, features = c("CD3D", "CD8A", "ITGAM", "CD14", "CD19", "NCAM1", "EOMES", "IL7R", "KLRB1", "ID2", "IL2RA", "THY1", "TBX21", "GATA3", "IL23R"), min.cutoff = "q9")
dev.off()


## Subset to cells with possible ILC identity (large bunch of clusters with scattered ILC marker expression, general lymphoid identity markers)
seurat.12.inter <- subset(seurat.12, idents = c(0, 1, 2, 4, 6, 7, 8, 9, 14))

## Cluster the cells
DefaultAssay(seurat.12.inter) <- "integrated"
seurat.12.inter <- FindNeighbors(seurat.12.inter, dims = 1:10)
seurat.12.inter <- FindClusters(seurat.12.inter, resolution = 0.5)

### Look at cluster IDs of the first 5 cells
head(Idents(seurat.12.inter), 5)

## Run non-linear dimensional reduction (UMAP/tSNE)
seurat.12.inter <- RunUMAP(seurat.12.inter, dims = 1:10)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12inter-initialUMAP.png")
DimPlot(seurat.12.inter, reduction = "umap", label = TRUE)
dev.off()
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12inter-initialUMAPcond.png")
DimPlot(seurat.12.inter, reduction = "umap", group.by = "condition")
dev.off()

## Find markers that differentiate clusters from one another (to help assign identities)
DefaultAssay(seurat.12.inter) <- "RNA"
seurat.12.inter.markers <- FindAllMarkers(seurat.12.inter, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.12.ilc <- seurat.12.inter.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.table(manual.curate.12.ilc, file="manual.curate.12.ilc.txt", row.names=TRUE, col.names=TRUE)

## Visualize specific markers to identify cell types
### ILC1s: CD127/IL7R+, CD161/KLRB1+, ID2+, T-bet/TBX21+, CD294/PTGDR2-, CD117/KIT-
## ILC2s: CD127/IL7R+, CD161/KLRB1+, ID2+, CD294/PTGDR2+, GATA3+, CD117/KIT+-
### ILC3s: CD127/IL7R+, CD161/KLRB1+, ID2+, CD294/PTGDR-, CD117/KIT+, IL23R+
### NKs: CD56/NCAM1+, EOMES+
### Ts: CD3/CD3D+ (definite), TBX21+
### Monocytes: CD14+
### Neutrophils, Basophils, Eosinophils, DCs: CD11b+
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12inter-initialmarkers.png", width= 700, height=480)
FeaturePlot(seurat.12.inter, features = c("CD3D", "CD8A", "ITGAM", "CD14", "CD19", "NCAM1", "EOMES", "IL7R", "KLRB1", "ID2", "IL2RA", "THY1", "TBX21", "GATA3", "IL23R"), min.cutoff = "q9")
dev.off()

## Assign identity to clusters
new.cluster.ids <- c("T cell", "T cell", "T cell", "NK cell or ILC", "T cell", "NK cell or ILC", "T cell", "Indeterminate", "T or NK cell", "Indeterminate")
names(new.cluster.ids) <- levels(seurat.12.inter)
seurat.12.inter <- RenameIdents(seurat.12.inter, new.cluster.ids)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12inter-labeledUMAP.png")
DimPlot(seurat.12.inter, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

## Subset to cells with possible ILC identity (NK or ILC, NK or T)
seurat.12.ilc <- subset(seurat.12.inter, idents = c("NK cell or ILC", "T or NK cell"))

## Cluster the cells
DefaultAssay(seurat.12.ilc) <- "integrated"
seurat.12.ilc <- FindNeighbors(seurat.12.ilc, dims = 1:10)
seurat.12.ilc <- FindClusters(seurat.12.ilc, resolution = 0.5)

### Look at cluster IDs of the first 5 cells
head(Idents(seurat.12.ilc), 5)

## Run non-linear dimensional reduction (UMAP/tSNE)
seurat.12.ilc <- RunUMAP(seurat.12.ilc, dims = 1:10)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12ilc-initialUMAP.png")
DimPlot(seurat.12.ilc, reduction = "umap", label = TRUE)
dev.off()
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12ilc-initialUMAPcond.png")
DimPlot(seurat.12.ilc, reduction = "umap", group.by = "condition")
dev.off()

## Find markers that differentiate clusters from one another (to help assign identities)
DefaultAssay(seurat.12.ilc) <- "RNA"
seurat.12.ilc.markers <- FindAllMarkers(seurat.12.ilc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.12.ilc <- seurat.12.ilc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)
write.table(manual.curate.12.ilc, file="manual.curate.12.ilc.txt", row.names=TRUE, col.names=TRUE)

## Visualize specific markers to identify cell types
### ILC1s: CD127/IL7R+, CD161/KLRB1+, ID2+, T-bet/TBX21+, CD294/PTGDR2-, CD117/KIT-
## ILC2s: CD127/IL7R+, CD161/KLRB1+, ID2+, CD294/PTGDR2+, GATA3+, CD117/KIT+-
### ILC3s: CD127/IL7R+, CD161/KLRB1+, ID2+, CD294/PTGDR-, CD117/KIT+, IL23R+
### NKs: CD56/NCAM1+, EOMES+
### Ts: CD3/CD3D+ (definite), TBX21+
### Monocytes: CD14+
### Neutrophils, Basophils, Eosinophils, DCs: CD11b+
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12ilc-initialmarkers.png", width= 700, height=480)
FeaturePlot(seurat.12.ilc, features = c("CD3D", "CD8A", "ITGAM", "CD14", "CD19", "NCAM1", "EOMES", "IL7R", "KLRB1", "ID2", "IL2RA", "THY1", "TBX21", "GATA3", "IL23R"), min.cutoff = "q9")
dev.off()

## Assign identity to clusters
## Cells do not have discernable levels of specific identifying markers; UMAP/ clustering did not separate T cells from NKs and ILCs by marker expression
## Manual marker ID
# new.cluster.ids <- c("")
# names(new.cluster.ids) <- levels(seurat.12.ilc)
# seurat.12.ilc <- RenameIdents(seurat.12.ilc, new.cluster.ids)
# png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12ilc-labeledUMAP.png")
# DimPlot(seurat.12.ilc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# dev.off()

## Identify conserved markers
seurat.12.interc.markers <- FindConservedMarkers(seurat.12.inter, ident.1 = c("NK cell or ILC", "T or NK cell"), grouping.var = "condition", verbose = FALSE)
write.table(seurat.12.interc.markers, file="conserved.12.inter.txt", row.names=TRUE, col.names=TRUE)

## Visualize specific immune markers by condition and cell type
markers.to.plot <- c("CD3D", "CD8A", "ITGAM", "CD14", "CD19", "NCAM1", "EOMES", "IL7R", "KLRB1", "ID2", "IL2RA", "THY1", "TBX21", "GATA3", "IL23R")
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat*12inter-immdotplot.png")
DotPlot(seurat.12.inter, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "condition") +
  RotatedAxis()
dev.off()

## Visualize specific microbiome interaction genes by condition and cell type
markers.to.plot <- c("NOD2", "CARD9", "ATG16L1*", "IRGM")
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12inter-microdotplot.png")
DotPlot(seurat.12.inter, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "condition") +
  RotatedAxis()
dev.off()

seurat.12.inter$celltype <- Idents(seurat.12.inter)
plots <- VlnPlot(seurat.12.inter, features = c("NOD2", "CARD9", "ATG16L1", "IRGM"), split.by = "condition", group.by = "celltype",
    pt.size = 0, combine = FALSE)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12inter-vinmicrogenes.png")
wrap_plots(plots = plots, ncol = 1)
dev.off()

## Visualize specific microbiome interaction genes in cowplot format by condition and cell type
## All cells
theme_set(theme_cowplot())
Idents(seurat.12) <- "condition"
avg.seurat.12 <- as.data.frame(log1p(AverageExpression(seurat.12, verbose = FALSE)$RNA))
avg.seurat.12$gene <- rownames(avg.seurat.12)

genes.to.label = c("NOD2", "CARD9", "ATG16L1", "IRGM") 
p9 <- ggplot(avg.seurat.12, aes(healthy, UC)) + geom_point() + ggtitle("All CD45+ cells")
p9 <- LabelPoints(plot = p9, points = genes.to.label, repel = TRUE)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12-diffmicrogenes.png")
p9
dev.off()

## Subsetted cells
Idents(seurat.12.ilc) <- "condition"
avg.ilcs <- as.data.frame(log1p(AverageExpression(seurat.12.ilc, verbose = FALSE)$RNA))
avg.ilcs$gene <- rownames(avg.ilcs)

genes.to.label = c("NOD2", "CARD9", "ATG16L1", "IRGM") 
p9 <- ggplot(avg.ilcs, aes(healthy, UC)) + geom_point() + ggtitle("Selected immune cells")
p9 <- LabelPoints(plot = p9, points = genes.to.label, repel = TRUE)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12ilc-diffmicrogenes.png")
p9
dev.off()



## Save objects for future recall
saveRDS(seurat.12, file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat12-clustered.rds")
saveRDS(seurat.12.inter, file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat12inter-clustered.rds")
saveRDS(seurat.12.ilc, file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat12ilc-clustered.rds")

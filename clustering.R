
# GSE 150050, colon tissues
## Cluster the cells
seurat.15 <- FindNeighbors(seurat.15, dims = 1:10)
seurat.15 <- FindClusters(seurat.15, resolution = 0.5)

### Look at cluster IDs of the first 5 cells
head(Idents(seurat.15), 5)

## Run non-linear dimensional reduction (UMAP/tSNE)
seurat.15 <- RunUMAP(seurat.15, dims = 1:10)
DimPlot(seurat.15, reduction = "umap", label = TRUE)

## Find markers that differentiate clusters from one another (to help assign identities)
seurat.15.markers <- FindAllMarkers(seurat.15, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.15 <- seurat.15.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)
write.table(manual.curate.15, file="manual.curate.15.txt", row.names=TRUE, col.names=TRUE)

## Visualize specific markers to identify immune cells 
### ILC1s: CD127/IL7R+, CD161/KLRB1+, CD117/KIT-, CD294/PTGDR2-, CD3-
### ILC2s: CD127/IL7R+, CD161/KLRB1+, CD294/PTGDR2+, GATA3+, CD3-
### ILC3s: CD127/IL7R+, CD161/KLRB1+, CD117/KIT+, CD3-
### NKs: CD56/NCAM1+,CD3-, EOMES+
### Ts: CD3/CD3D+ (definite)
FeaturePlot(seurat.15, features = c("CD3D", "NCAM1", "EOMES", "IL7R", "PTGDR2", "KIT", "KLRB1", "GATA3"), min.cutoff = "q9")

## Assign cell type identity to clusters
new.cluster.ids <- c("ILC3", "T cell", "ILC1", "ILC3", "Indeterminate", "Indeterminate", "Indeterminate ILC", "T cell", "T cell", "NK cell", "Indeterminate", "ILC2", "Indeterminate", "Indeterminate ILC")
names(new.cluster.ids) <- levels(seurat.15)
seurat.15 <- RenameIdents(seurat.15, new.cluster.ids)
DimPlot(seurat.15, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## Visualize specific markers by cluster
Idents(seurat.15) <- factor(Idents(seurat.15), levels = c("ILC3", "T cell", "ILC1", "ILC3", "Indeterminate", "Indeterminate", "Indeterminate ILC", "T cell", "T cell", "NK cell", "Indeterminate", "ILC2", "Indeterminate", "Indeterminate ILC"))
markers.to.plot <- c("CD3D", "NCAM1", "EOMES", "IL7R", "PTGDR2", "KIT", "KLRB1", "GATA3")
DotPlot(seurat.15, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()

## Subset to cells with determinate and non-T cell identities
seurat.15.ilc <- subset(seurat.15, idents = c("ILC1", "ILC2", "ILC3", "NK cell", "Indeterminate ILC"))

## Cluster subset
seurat.15.ilc <- FindNeighbors(seurat.15.ilc, dims = 1:10)
seurat.15.ilc <- FindClusters(seurat.15.ilc, resolution = 0.5)

## Run non-linear dimensional reduction (UMAP/tSNE)
seurat.15.ilc <- RunUMAP(seurat.15.ilc, dims = 1:10)
DimPlot(seurat.15.ilc, reduction = "umap", label = TRUE)

## Find markers that differentiate clusters from one another (to help assign identities)
seurat.15.ilc.markers <- FindAllMarkers(seurat.15.ilc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.15.ilc <- seurat.15.ilc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)
write.table(manual.curate.15.ilc, file="manual.curate.15.ilc.txt", row.names=TRUE, col.names=TRUE)

## Visualize specific markers to identify immune cell types
### ILC1s: CD127/IL7R+, CD161/KLRB1+, CD117/KIT-, CD294/PTGDR2-
### ILC2s: CD127/IL7R+, CD161/KLRB1+, CD117/KIT+-, CD294/PTGDR2+
### ILC3s: CD127/IL7R+, CD161/KLRB1+, CD117/KIT+, CD294/PTGDR-
### NKs: CD56/NCAM1+, EOMES+
FeaturePlot(seurat.15.ilc, features = c("NCAM1", "EOMES", "IL7R", "PTGDR2", "KIT", "KLRB1", "GATA3"), min.cutoff = "q9")

## Assign cell type identity to clusters
new.cluster.ids <- c("ILC3", "ILC3", "ILC1", "ILC3", "NK", "NK", "ILC2", "ILC2", "Indeterminate")
names(new.cluster.ids) <- levels(seurat.15.ilc)
seurat.15.ilc <- RenameIdents(seurat.15.ilc, new.cluster.ids)
DimPlot(seurat.15.ilc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## Visualize specific markers by cell type
Idents(seurat.15.ilc) <- factor(Idents(seurat.15.ilc), levels = c("ILC3", "ILC3", "ILC1", "ILC3", "NK", "NK", "ILC2", "ILC2", "Indeterminate"))
markers.to.plot <- c("NCAM1", "EOMES", "IL7R", "PTGDR2", "KIT", "KLRB1", "GATA3")
DotPlot(seurat.15.ilc, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()

## Save objects for future recall
saveRDS(seurat.15, file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat15-clustered.rds")
saveRDS(seurat.15.ilc, file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat15ilc-clustered.rds")


# GSE 185224, normal gut
## Cluster the cells
seurat.18 <- FindNeighbors(seurat.18, dims = 1:10)
seurat.18 <- FindClusters(seurat.18, resolution = 0.5)

### Look at cluster IDs of the first 5 cells
head(Idents(seurat.18), 5)

## Run non-linear dimensional reduction (UMAP/tSNE)
seurat.18 <- RunUMAP(seurat.18, dims = 1:10)
DimPlot(seurat.18, reduction = "umap", label = TRUE)

## Find markers that differentiate clusters from one another (to help assign identities)
seurat.18.markers <- FindAllMarkers(seurat.18, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.18 <- seurat.18.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.table(manual.curate.18, file="manual.curate.18.txt", row.names=TRUE, col.names=TRUE)

## Assign cell type identity to clusters
new.cluster.ids <- c("Distal enterocyte", "Indeterminate enterocyte", "Undifferentiated", "Proximal enterocyte", "Proximal enterocyte", "Paneth", "Undifferentiated", "Intestinal globlet", "Undifferentiated", "Distal enterocyte", "Immune", "Proximal enterocyte", "Indeterminate enterocyte", "Paneth", "Immune", "Paneth", "Enteroendocrine")
names(new.cluster.ids) <- levels(seurat.18)
seurat.18 <- RenameIdents(seurat.18, new.cluster.ids)
DimPlot(seurat.18, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## Subset clusters with immune identity
seurat.18.immune <- subset(seurat.18, idents = "Immune")

## Cluster immune cell subset
seurat.18.immune <- FindNeighbors(seurat.18.immune, dims = 1:10)
seurat.18.immune <- FindClusters(seurat.18.immune, resolution = 0.5)

### Look at cluster IDs of the first 5 cells
head(Idents(seurat.18.immune), 5)

## Run non-linear dimensional reduction (UMAP/tSNE)
seurat.18.immune <- RunUMAP(seurat.18.immune, dims = 1:10)
DimPlot(seurat.18.immune, reduction = "umap", label = TRUE)

## Find markers that differentiate clusters from one another (to help assign identities)
seurat.18.immune.markers <- FindAllMarkers(seurat.18.immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.18.immune <- seurat.18.immune.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.table(manual.curate.18.immune, file="manual.curate.18.immune.txt", row.names=TRUE, col.names=TRUE)

## Visualize specific markers to identify immune cells from a predominant epithelial population
## ILCs (CD127/IL7R+, CD117/KIT +, CD161/KLRB1+); no distinct signatures found
FeaturePlot(seurat.18.immune, features = c("CD4", "CD8A", "IL7R", "PTGDR2", "KIT", "CD14", "CD19", "KLRB1"), min.cutoff = "q9")

new.cluster.ids <- c("")
names(new.cluster.ids) <- levels(seurat.18.immune)
seurat.18 <- RenameIdents(seurat.18, new.cluster.ids)
DimPlot(seurat.18.immune, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## Plot specific markers
Idents(seurat.18.immune) <- factor(Idents(seurat.18.immune), levels = c(""))
markers.to.plot <- c("CD4", "CD8A", "IL7R", "PTGDR2", "KIT", "CD14", "CD19", "KLRB1")
DotPlot(seurat.18.immune, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()

# Save objects for future recall
saveRDS(seurat.18, file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat18-clustered.rds")
saveRDS(seurat.18.immune, file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat18immune-clustered.rds")


# GSE 125527, IBD gut
# Intestinal immune cells
## Cluster the cells
seurat.12.int <- FindNeighbors(seurat.12.int, dims = 1:10)
seurat.12.int <- FindClusters(seurat.12.int, resolution = 0.5)

### Look at cluster IDs of the first 5 cells
head(Idents(seurat.12.int), 5)

## Run non-linear dimensional reduction (UMAP/tSNE)
seurat.12.int <- RunUMAP(seurat.12.int, dims = 1:10)
DimPlot(seurat.12.int, reduction = "umap", label = TRUE)

## Find markers that differentiate clusters from one another (to help assign identities)
seurat.12.int.markers <- FindAllMarkers(seurat.12.int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.12.int <- seurat.12.int.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)
write.table(manual.curate.12.int, file="manual.curate.12.int.txt", row.names=TRUE, col.names=TRUE)

## Visualize specific markers to identify immune cells 
### ILC1s: CD127/IL7R+, CD161/KLRB1+, CD117/KIT-, CD294/PTGDR2-, CD3-
### ILC2s: CD127/IL7R+, CD161/KLRB1+, CD294/PTGDR2+, GATA3+, CD3-
### ILC3s: CD127/IL7R+, CD161/KLRB1+, CD117/KIT+, CD3-
### NKs: CD56/NCAM1+,CD3-, EOMES+
### Ts: CD3/CD3D+ (definite)
### Bs: 
### Myeloids: 
FeaturePlot(seurat.12.int, features = c("CD3D", "NCAM1", "EOMES", "IL7R", "PTGDR2", "KIT", "KLRB1", "GATA3"), min.cutoff = "q9")

## Assign cell type identity to clusters
new.cluster.ids <- c("")
names(new.cluster.ids) <- levels(seurat.12.int)
seurat.12.int <- RenameIdents(seurat.12.int, new.cluster.ids)
DimPlot(seurat.12.int, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## Visualize specific markers by cluster
Idents(seurat.12.int) <- factor(Idents(seurat.12.int), levels = c(""))
markers.to.plot <- c("")
DotPlot(seurat.12.int, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()

## Subset to cells with determinate and non-T cell identities
seurat.12.int.ilc <- subset(seurat.12.int, idents = c(""))

## Cluster subset
seurat.12.int.ilc <- FindNeighbors(seurat.12.int.ilc, dims = 1:10)
seurat.12.int.ilc <- FindClusters(seurat.12.int.ilc, resolution = 0.5)

## Run non-linear dimensional reduction (UMAP/tSNE)
seurat.12.int.ilc <- RunUMAP(seurat.12.int.ilc, dims = 1:10)
DimPlot(seurat.12.int.ilc, reduction = "umap", label = TRUE)

## Find markers that differentiate clusters from one another (to help assign identities)
seurat.12.int.ilc.markers <- FindAllMarkers(seurat.12.int.ilc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.12.int.ilc <- seurat.12.int.ilc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 3, order_by = avg_log2FC)
write.table(manual.curate.12.int.ilc, file="manual.curate.12.int.ilc.txt", row.names=TRUE, col.names=TRUE)

## Visualize specific markers to identify immune cell types
### ILC1s: CD127/IL7R+, CD161/KLRB1+, CD117/KIT-, CD294/PTGDR2-
### ILC2s: CD127/IL7R+, CD161/KLRB1+, CD117/KIT+-, CD294/PTGDR2+
### ILC3s: CD127/IL7R+, CD161/KLRB1+, CD117/KIT+, CD294/PTGDR-
### NKs: CD56/NCAM1+, EOMES+
FeaturePlot(seurat.12.int.ilc, features = c("NCAM1", "EOMES", "IL7R", "PTGDR2", "KIT", "KLRB1", "GATA3"), min.cutoff = "q9")

## Assign cell type identity to clusters
new.cluster.ids <- c("")
names(new.cluster.ids) <- levels(seurat.12.int.ilc)
seurat.12.int.ilc <- RenameIdents(seurat.12.int.ilc, new.cluster.ids)
DimPlot(seurat.12.int.ilc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## Visualize specific markers by cell type
Idents(seurat.12.int.ilc) <- factor(Idents(seurat.12.int.ilc), levels = c(""))
markers.to.plot <- c("NCAM1", "EOMES", "IL7R", "PTGDR2", "KIT", "KLRB1", "GATA3")
DotPlot(seurat.12.int.ilc, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()

## Save objects for future recall
saveRDS(seurat.12.int, file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat12int-clustered.rds")
saveRDS(seurat.12.int.ilc, file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/rds/seurat12intilc-clustered.rds")


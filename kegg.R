# GSE 150050
## Extract differentially expressed features matrices
## All-cluster matrix
im.15.markers <- FindAllMarkers(seurat.15, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.15 <- seurat.15.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
im.15.markers <- im.15.markers %>% select(p_val, avg_log2FC)
im.15.markers <- rownames_to_column(im.15.markers, "gene_id")

## Individual cluster matrices
im.cluster0.markers <- FindMarkers(seurat.15.int, ident.1 = 2, min.pct = 0.25)
im.cluster0.markers <- im.cluster0.markers %>% select(p_val, avg_log2FC)
im.cluster0.markers <- rownames_to_column(im.cluster0.markers, "gene_id")

### ILC versus other clusters matrix
im.ilc.markers <- FindMarkers(seurat.12.int.ilc, ident.1 = c("ILC1", "ILC2", "ILC3"), min.pct = 0.25)
im.ilc.markers <- im.ilc.markers %>% select(p_val, avg_log2FC)
im.ilc.markers <- rownames_to_column(im.ilc.markers, "gene_id")

### ILC subsets versus each other matrices
im.ilc1.markers <- FindMarkers(seurat.12.int.ilc, ident.1 = "ILC1", min.pct = 0.25)
im.ilc1.markers <- im.ilc1.markers %>% select(p_val, avg_log2FC)
im.ilc1.markers <- rownames_to_column(im.ilc1.markers, "gene_id")

im.ilc2.markers <- FindMarkers(seurat.12.int.ilc, ident.1 = "ILC2", min.pct = 0.25)
im.ilc2.markers <- im.ilc2.markers %>% select(p_val, avg_log2FC)
im.ilc2.markers <- rownames_to_column(im.ilc2.markers, "gene_id")

im.ilc3.markers <- FindMarkers(seurat.12.int.ilc, ident.1 = "ILC3", min.pct = 0.25)
im.ilc3.markers <- im.ilc3.markers %>% select(p_val, avg_log2FC)
im.ilc3.markers <- rownames_to_column(im.ilc3.markers, "gene_id")

# GSE 125527 
## Intestinal immune cells
## Extract differentially expressed features matrices
## All-cluster matrix
int.12.markers <- FindAllMarkers(seurat.12.int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.12.int <- seurat.12.int.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
int.12.markers <- int.12.markers %>% select(p_val, avg_log2FC)
int.12.markers <- rownames_to_column(int.12.markers, "gene_id")

## Individual cluster matrices
int.cluster0.markers <- FindMarkers(seurat.12.int.ilc, ident.1 = 2, min.pct = 0.25)
int.cluster0.markers <- int.cluster0.markers %>% select(p_val, avg_log2FC)
int.cluster0.markers <- rownames_to_column(int.cluster0.markers, "gene_id")

### ILC versus other clusters matrix
int.ilc.markers <- FindMarkers(seurat.12.int.ilc, ident.1 = c("ILC1", "ILC2", "ILC3"), min.pct = 0.25)
int.ilc.markers <- int.ilc.markers %>% select(p_val, avg_log2FC)
int.ilc.markers <- rownames_to_column(int.ilc.markers, "gene_id")

### ILC subsets versus each other matrices
int.ilc1.markers <- FindMarkers(seurat.12.int.ilc, ident.1 = "ILC1", min.pct = 0.25)
int.ilc1.markers <- int.ilc1.markers %>% select(p_val, avg_log2FC)
int.ilc1.markers <- rownames_to_column(int.ilc1.markers, "gene_id")

int.ilc2.markers <- FindMarkers(seurat.12.int.ilc, ident.1 = "ILC2", min.pct = 0.25)
int.ilc2.markers <- int.ilc2.markers %>% select(p_val, avg_log2FC)
int.ilc2.markers <- rownames_to_column(int.ilc2.markers, "gene_id")

int.ilc3.markers <- FindMarkers(seurat.12.int.ilc, ident.1 = "ILC3", min.pct = 0.25)
int.ilc3.markers <- int.ilc3.markers %>% select(p_val, avg_log2FC)
int.ilc3.markers <- rownames_to_column(int.ilc3.markers, "gene_id")


# PBMCs
## Extract differentially expressed features matrices
## All-cluster matrix
pbmc.12.markers <- FindAllMarkers(seurat.12, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.12 <- seurat.12.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
pbmc.12.markers <- pbmc.12.markers %>% select(p_val, avg_log2FC)
pbmc.12.markers <- rownames_to_column(pbmc.12.markers, "gene_id")

## Individual cluster 
cluster0.markers <- FindMarkers(seurat.12, ident.1 = 2, min.pct = 0.25)
cluster0.markers <- cluster0.markers %>% select(p_val, avg_log2FC)
cluster0.markers <- rownames_to_column(cluster0.markers, "gene_id")

### ILC versus other clusters 
ilc.markers <- FindMarkers(seurat.12, ident.1 = c("ILC1", "ILC2", "ILC3"), min.pct = 0.25)
ilc.markers <- ilc.markers %>% select(p_val, avg_log2FC)
ilc.markers <- rownames_to_column(ilc.markers, "gene_id")

### ILC subsets 
ilc1.markers <- FindMarkers(seurat.12, ident.1 = "ILC1", min.pct = 0.25)
ilc1.markers <- ilc1.markers %>% select(p_val, avg_log2FC)
ilc1.markers <- rownames_to_column(ilc1.markers, "gene_id")

ilc2.markers <- FindMarkers(seurat.12, ident.1 = "ILC2", min.pct = 0.25)
ilc2.markers <- ilc2.markers %>% select(p_val, avg_log2FC)
ilc2.markers <- rownames_to_column(ilc2.markers, "gene_id")

ilc3.markers <- FindMarkers(seurat.12, ident.1 = "ILC3", min.pct = 0.25)
ilc3.markers <- ilc3.markers %>% select(p_val, avg_log2FC)
ilc3.markers <- rownames_to_column(ilc3.markers, "gene_id")

## Differences by condition
### All cluster matrix, by condition
hu.12.markers <- FindMarkers(seurat.12, ident.1 = "healthy", ident.2 = "UC", verbose = FALSE)
hu.12.markers <- hu.12.markers %>% select(p_val, avg_log2FC)
hu.12.markers <- rownames_to_column(hu.12.markers, "gene_id")

### ILC matrix by condition
hu.ilc.markers <- FindMarkers(seurat.12.ilc, ident.1 = "healty", ident.2 = "UC",  min.pct = 0.25)
hu.ilc.markers <- hu.ilc.markers %>% select(p_val, avg_log2FC)
hu.ilc.markers <- rownames_to_column(hu.ilc.markers, "gene_id")

### ILC subsets versus the same subset by condition 


# Run pathfindR
## GSE 150050
im.15.path <- run_pathfindR(im.15.markers)
im.cluster0.path <- run_pathfindR(im.cluster0.markers)
im.ilc.path <- run_pathfindR(im.ilc.markers)
im.ilc1.path <- run_pathfindR(im.ilc1.markers)
im.ilc2.path <- run_pathfindR(im.ilc2.markers)
im.ilc3.path <- run_pathfindR(im.ilc3.markers)

## GSE 125527
## Intestinal immune cells
int.12.path <- run_pathfindR(int.12.markers)
int.cluster0.path <- run_pathfindR(int.cluster0.markers)
int.ilc.path <- run_pathfindR(int.ilc.markers)
int.ilc1.path <- run_pathfindR(int.ilc1.markers)
int.ilc2.path <- run_pathfindR(int.ilc2.markers)
int.ilc3.path <- run_pathfindR(int.ilc3.markers)

## PBMCs 
pbmc.12.path <- run_pathfindR(pbmc.12.markers)
cluster0.path <- run_pathfindR(cluster0.markers)
ilc.path <- run_pathfindR(ilc.markers)
ilc1.path <- run_pathfindR(ilc1.markers)
ilc2.path <- run_pathfindR(ilc2.markers)
ilc3.path <- run_pathfindR(ilc3.markers)

## PBMCs by condition
hu.12.path <- run_pathfindR(hu.12.markers)
hu.ilc.path <- run_pathfindR(hu.ilc.markers)


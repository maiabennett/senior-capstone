# GSE 150050
## All-cluster matrix
im.15.markers <- FindAllMarkers(seurat.15, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.15 <- im.15.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
im.15.markers <- im.15.markers %>% select(p_val, avg_log2FC)
im.15.markers <- rownames_to_column(im.15.markers, "gene_id")
# Replace 0s
im.15.markers$p_val[im.15.markers$p_val=="0"]<-1.0e-302
im.15.markers <- im.15.markers[c("gene_id", "avg_log2FC", "p_val")]

### ILC versus other clusters matrix
im.ilc.markers <- FindMarkers(seurat.15.ilc, ident.1 = c("ILC1", "ILC2", "ILC3"), min.pct = 0.25)
im.ilc.markers <- im.ilc.markers %>% select(p_val, avg_log2FC)
im.ilc.markers <- rownames_to_column(im.ilc.markers, "gene_id")

### ILC subsets versus each other matrices
im.ilc1.markers <- FindMarkers(seurat.15.ilc, ident.1 = "ILC1", indent.2 = c("ILC2", "ILC3"), min.pct = 0.25)
im.ilc1.markers <- im.ilc1.markers %>% select(p_val, avg_log2FC)
im.ilc1.markers <- rownames_to_column(im.ilc1.markers, "gene_id")
# Replace 0s
im.ilc1.markers$p_val[im.ilc1.markers$p_val=="0"]<-1.0e-302

im.ilc2.markers <- FindMarkers(seurat.15.ilc, ident.1 = "ILC2", indent.2 = c("ILC1", "ILC3"), min.pct = 0.25)
im.ilc2.markers <- im.ilc2.markers %>% select(p_val, avg_log2FC)
im.ilc2.markers <- rownames_to_column(im.ilc2.markers, "gene_id")
# Replace 0s
im.ilc2.markers$p_val[im.ilc2.markers$p_val=="0"]<-1.0e-302

im.ilc3.markers <- FindMarkers(seurat.15.ilc, ident.1 = "ILC3", indent.2 = c("ILC1", "ILC2"), min.pct = 0.25)
im.ilc3.markers <- im.ilc3.markers %>% select(p_val, avg_log2FC)
im.ilc3.markers <- rownames_to_column(im.ilc3.markers, "gene_id")
# Replace 0s
im.ilc3.markers$p_val[im.ilc3.markers$p_val=="0"]<-1.0e-302

# GSE 125527 
## Intestinal immune cells
## All-cluster matrix
int.12.markers <- FindAllMarkers(seurat.12.int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.12.int <- seurat.12.int.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
int.12.markers <- int.12.markers %>% select(p_val, avg_log2FC)
int.12.markers <- rownames_to_column(int.12.markers, "gene_id")
# Replace 0s
int.12.markers$p_val[int.12.markers$p_val=="0"]<-1.0e-302

### ILC versus other clusters matrix
int.ilc.markers <- FindMarkers(seurat.12.int.ilc, ident.1 = "Contains ILCs", min.pct = 0.25)
int.ilc.markers <- int.ilc.markers %>% select(p_val, avg_log2FC)
int.ilc.markers <- rownames_to_column(int.ilc.markers, "gene_id")
# Replace 0s
int.ilc.markers$p_val[int.ilc.markers$p_val=="0"]<-1.0e-302


# PBMCs
## All-cluster matrix
pbmc.12.markers <- FindAllMarkers(seurat.12, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
manual.curate.12 <- seurat.12.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)
pbmc.12.markers <- pbmc.12.markers %>% select(p_val, avg_log2FC)
pbmc.12.markers <- rownames_to_column(pbmc.12.markers, "gene_id")

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


# Run pathfindR
## GSE 150050
im.15.path <- run_pathfindR(im.15.markers, output_dir = "~/pathfindR/seurat15_all")
im.ilc.path <- run_pathfindR(im.ilc.markers, output_dir = "~/pathfindR/seurat15_ilc")
im.ilc1.path <- run_pathfindR(im.ilc1.markers, output_dir = "~/pathfindR/seurat15_ilc1")
im.ilc2.path <- run_pathfindR(im.ilc2.markers, output_dir = "~/pathfindR/seurat15_ilc2")
im.ilc3.path <- run_pathfindR(im.ilc3.markers, output_dir = "~/pathfindR/seurat15_ilc2")

## GSE 125527
## Intestinal immune cells
int.12.path <- run_pathfindR(int.12.markers, output_dir = "~/pathfindR/seurat12int_all")
int.ilc.path <- run_pathfindR(int.ilc.markers, output_dir = "~/pathfindR/seurat12int_ilc")

## PBMCs 
pbmc.12.path <- run_pathfindR(pbmc.12.markers)
ilc.path <- run_pathfindR(ilc.markers)
ilc1.path <- run_pathfindR(ilc1.markers)
ilc2.path <- run_pathfindR(ilc2.markers)
ilc3.path <- run_pathfindR(ilc3.markers)

## PBMCs by condition
hu.12.path <- run_pathfindR(hu.12.markers)
hu.ilc.path <- run_pathfindR(hu.ilc.markers)


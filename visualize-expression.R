## Visualize specific microbiome interaction genes by condition and cell type
## GSE 150050
markers.to.plot <- c("NOD2", "CARD9", "ATG16L1*", "IRGM")
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat15-microdotplot.png")
DotPlot(seurat.15.ilc, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()
dev.off()

seurat.15.ilc$celltype <- Idents(seurat.15.ilc)
plots <- VlnPlot(seurat.15.ilc, features = c("NOD2", "CARD9", "ATG16L1", "IRGM"), group.by = "celltype",
    pt.size = 0, combine = FALSE, y.max = 2)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat15ilc-vinmicrogenes.png")
wrap_plots(plots = plots, ncol = 1)
dev.off()


## GSE 125527 intestines
markers.to.plot <- c("NOD2", "CARD9", "ATG16L1*", "IRGM")
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12intilc-microdotplot.png")
DotPlot(seurat.12.int.ilc, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8) +
  RotatedAxis()
dev.off()

seurat.12.int.ilc$celltype <- Idents(seurat.12.int.ilc)
plots <- VlnPlot(seurat.12.int.ilc, features = c("NOD2", "CARD9", "ATG16L1", "IRGM"), group.by = "celltype",
    pt.size = 0, combine = FALSE, y.max = 2)
png(file = "C:/Users/Me/OneDrive - University of Nebraska at Omaha/Administrative/Documents/Senior Project/senior-capstone/graphics/seurat12intilc-vinmicrogenes.png")
wrap_plots(plots = plots, ncol = 1)
dev.off()

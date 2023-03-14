
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
pbmc.12.path <- run_pathfindR(pbmc.12.markers, output_dir = "~/pathfindR/seurat12_all")
inter.ilc.path <- run_pathfindR(inter.ilc.markers, output_dir = "~/pathfindR/seurat12_ilc")

## PBMCs by condition
hu.12.path <- run_pathfindR(hu.12.markers, output_dir = "~/pathfindR/seurat12_conditions")
hu.ilc.path <- run_pathfindR(hu.ilc.markers, output_dir = "~/pathfindR/seurat12_ilc_conditions")


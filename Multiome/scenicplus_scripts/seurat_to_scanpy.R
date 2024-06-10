library(Seurat)
library(SeuratDisk)

seurat_obj <- readRDS("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240517_seuratv5_annotated.rds")
# Subset for 10,000 monocytes
# cell_barcodes <- read.table('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_10000_cells_20_topics/cell_barcodes.txt')$V1
cell_barcodes <- read.table("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/matrices/monocyte/barcodes.tsv.gz", sep='\t')

# Subest seurat object for cell barcodes
seurat_obj_subset <- seurat_obj[,colnames(seurat_obj) %in% cell_barcodes$V1]

# Write seurat to 10X format 
library(DropletUtils)
write10xCounts(x = seurat_obj_subset@assays$SCT@counts, path = '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/mono_all/seurat_10x_format')
library(Seurat)
library(SeuratDisk)

seurat_obj <- readRDS("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240223_seuratv5_normalized.rds")
cell_barcodes <- read.table('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_5000_cell_20_topics/cell_barcodes.txt')$V1

# Subest seurat object for cell barcodes
seurat_obj_subset <- seurat_obj[,colnames(seurat_obj) %in% cell_barcodes]

# THIS CODE IS NOT USED ANYMORE, 10X COUNT FORMAT IS ENOUGH
# # Resolves issue: https://github.com/mojaveazure/seurat-disk/issues/27
# seurat_obj_subset[["RNA"]] <- as(object = seurat_obj_subset[["RNA"]], Class = "Assay")
# seurat_obj_subset[["MJ"]] <- as(object = seurat_obj_subset[["MJ"]], Class = "Assay")

# SaveH5Seurat(seurat_obj_subset, filename = "/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/seurat_obj_scaled_data.h5Seurat", overwrite=TRUE)
# Convert("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/seurat_obj_scaled_data.h5Seurat", dest="h5ad")



# Write seurat to 10X format 
library(DropletUtils)
write10xCounts(x = seurat_obj_subset@assays$SCT@counts, path = '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/seurat_10x_format')
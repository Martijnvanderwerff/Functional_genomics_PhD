library(Seurat)
library(SeuratDisk)

seurat_obj <- readRDS("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240223_seuratv5_normalized.rds")

# Resolves issue: https://github.com/mojaveazure/seurat-disk/issues/27
seurat_obj[["RNA"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay")
seurat_obj[["MJ"]] <- as(object = seurat_obj[["MJ"]], Class = "Assay")

SaveH5Seurat(seurat_obj, filename = "/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/seurat_obj.h5Seurat", overwrite=TRUE)
seurat_to_scanpy <- Convert(seurat_obj, to = "anndata", filename = "/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/seurat_preprocessed.h5ad")
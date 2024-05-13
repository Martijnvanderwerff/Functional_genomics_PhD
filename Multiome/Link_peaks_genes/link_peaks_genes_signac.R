library(Seurat)
library(stringr)
library(Signac)
library(plyranges)
library(future)
plan("multicore", workers=4)

options(future.globals.maxSize = 200 * 1024 ^ 3) # for 50 Gb RAM

seurat_obj <- readRDS('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/seurat_preprocess_samples/objects/mo_all_20240223_seuratv5_normalized.rds')
# Get monocytes
seurat_obj <- seurat_obj[,seurat_obj@meta.data$celltype_imputed_lowerres == 'monocyte']

signac_peaks <- readRDS('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/cpeaks_peak_calling/signac/rounded/mo_cpeaks_filtered_percelltypemajor_wstatus_1_64.rds')
signac_peaks <- signac_peaks$monocyte

# 122,698 cells
common_barcodes <- intersect(names(signac_peaks$barcode), names(seurat_obj$barcode))

# Filter on common barcodes
seurat_obj_common <- seurat_obj[,colnames(seurat_obj) %in% common_barcodes]
signac_peaks_common <- signac_peaks[,colnames(signac_peaks) %in% common_barcodes]

# Extract assay from Seurat object of the monocytes
signac_peak_assay_data <- GetAssayData(signac_peaks_common, slot="data")
# Subset peaks for chr6
#chr6_peaks <- signac_peak_assay_data[startsWith(rownames(signac_peak_assay_data), "chr6"),]
# Create new Chromatin assay to add to seurat_obj_common
chromatinassay <- CreateChromatinAssay(signac_peak_assay_data)
# Add chromatin assay to existing seurat object with common cells
seurat_obj_common[['peaks']] <- chromatinassay

# Region stats
library(BSgenome.Hsapiens.UCSC.hg38)
seurat_obj_common[['peaks']] <- RegionStats(
  object = seurat_obj_common[['peaks']],
  genome = BSgenome.Hsapiens.UCSC.hg38
)

library(EnsDb.Hsapiens.v86)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
DefaultAssay(seurat_obj_common) <- 'peaks'
Annotation(seurat_obj_common) <- annotations

saveRDS(seurat_obj_common, '/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/objects/multimodal_seurat_obj_monocytes.rds')

# Look for chr6 genes
genes_chr6 <- read.csv('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/chr6_mono_genes.tsv', sep='\t', header=F)$V1

# Subset for donors
seurat_donor_objects <- SplitObject(seurat_obj_common, split.by="unconfined_best_match_sample")
saveRDS(seurat_donor_objects, "/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/objects/multimodal_seurat_obj_monocytes_per_donor.rds")
# Test first donor 
first_donor <- seurat_donor_objects[[1]]

seurat_donor_objects <- readRDS( "/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/objects/multimodal_seurat_obj_monocytes_per_donor.rds")

for (donor in names(seurat_donor_objects)) {
  saveRDS(seurat_donor_objects[[donor]], paste0('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/objects/multimodel_seurat_obj_all_donors_', donor, '.rds'))
}

for (donor in names(seurat_donor_objects)) {
  
  split_condition <- SplitObject(seurat_donor_objects[[donor]], split.by="")
  linked_peaks_ut <- LinkPeaks(split_condition[[]], peak.assay = 'peaks', expression.assay = 'SCT', genes.use=genes_chr6, method='spearman', distance=100000)
  linked_peaks_stim <- LinkPeaks(split_condition[[]], peak.assay = 'peaks', expression.assay = 'SCT', genes.use=genes_chr6, method='spearman', distance=100000)

  links_ut <- linked_peaks_ut@assays$peaks@links
  links_stim <- linked_peaks_stim@assays$peaks@links

  saveRDS(links, paste0("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/objects/region_gene_links_per_donor_ut_stim/links_monocytes_", donor, "_UT.rds"))
  saveRDS(links_stim, paste0("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/objects/region_gene_links_per_donor_ut_stim/links_monocytes_", donor, "_24hCA.rds"))
}

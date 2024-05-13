library(Seurat)
library(stringr)
library(Signac)
library(plyranges)
library(future)
plan("multicore", workers=4)

options(future.globals.maxSize = 200 * 1024 ^ 3) # for 50 Gb RAM

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args[1])
input_object = readRDS(paste0('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/objects/multimodel_seurat_obj_all_donors/multimodel_seurat_obj_all_donors_', args[1], '.rds'))
print(input_object)
genes_chr6 <- read.csv('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/chr6_mono_genes.tsv', sep='\t', header=F)$V1

if ("UT" %in% input_object@meta.data$confined_condition && "24hCa" %in% input_object@meta.data$confined_condition) {
    split_condition <- SplitObject(input_object, split.by="confined_condition")
    linked_peaks_ut <- LinkPeaks(split_condition[["UT"]], peak.assay = 'peaks', expression.assay = 'SCT', genes.use=genes_chr6, method='spearman', distance=100000)
    linked_peaks_stim <- LinkPeaks(split_condition[["24hCa"]], peak.assay = 'peaks', expression.assay = 'SCT', genes.use=genes_chr6, method='spearman', distance=100000)

    links_ut <- linked_peaks_ut@assays$peaks@links
    links_stim <- linked_peaks_stim@assays$peaks@links

    saveRDS(links_ut, paste0("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/objects/region_gene_links_per_donor_ut_stim/links_monocytes_", args[1], "_UT.rds"))
    saveRDS(links_stim, paste0("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/objects/region_gene_links_per_donor_ut_stim/links_monocytes_", args[1], "_24hCa.rds"))
}  else if (length(unique(input_object@meta.data$confined_condition)) == 1 && unique(input_object@meta.data$confined_condition) == "UT") {
    
    links <- LinkPeaks(input_object, peak.assay = 'peaks', expression.assay = 'SCT', genes.use=genes_chr6, method='spearman', distance=100000)
    links_ut <- links@assays$peaks@links
    saveRDS(links_ut, paste0("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/objects/region_gene_links_per_donor_ut_stim/links_monocytes_", args[1], "_UT.rds"))
} else {
    links <- LinkPeaks(input_object, peak.assay = 'peaks', expression.assay = 'SCT', genes.use=genes_chr6, method='spearman', distance=100000)
    links_stim <- links@assays$peaks@links
    saveRDS(links_stim, paste0("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/objects/region_gene_links_per_donor_ut_stim/links_monocytes_", args[1], "_24hCa.rds"))
}
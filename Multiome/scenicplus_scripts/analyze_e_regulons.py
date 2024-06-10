import mudata
import pandas as pd

# Load data
scplus_mdata = mudata.read('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/test_10000_monocytes_20_topics/outs/scplusmdata.h5mu')
directed_e_regulons = scplus_mdata.uns["direct_e_regulon_metadata"]
extended_e_regulons = scplus_mdata.uns["extended_e_regulon_metadata"]

directed_e_regulons.to_csv('direct_e_regulons.csv', index=None)
extended_e_regulons.to_csv('extended_e_regulons.csv', index=None)

from scenicplus.RSS import (regulon_specificity_scores, plot_rss)

rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scRNA_counts:second_most_contributing_topic",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC"]
)

import scanpy as sc
import anndata
eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,
)

eRegulon_gene_AUC.obs = scplus_mdata.obs.loc[eRegulon_gene_AUC.obs_names]

sc.pp.neighbors(eRegulon_gene_AUC, use_rep = "X")

sc.tl.umap(eRegulon_gene_AUC)

sc.pl.umap(eRegulon_gene_AUC, color = "scRNA_counts:celltype_imputed", save='_pca_enrichment_scores_col_subcelltype.pdf')

from scenicplus.plotting.dotplot import heatmap_dotplot

heatmap_dotplot(
    scplus_mudata = scplus_mdata,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scRNA_counts:condition_final",
    eRegulon_metadata_key = "direct_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    orientation = "horizontal",
    figsize = (20, 10),
    save='process_data_files/heatmap_enrichment_condition.pdf',
)

import scanpy as sc
import anndata
eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,
)

eRegulon_gene_AUC.obs = scplus_mdata.obs.loc[eRegulon_gene_AUC.obs_names]

sc.pp.neighbors(eRegulon_gene_AUC, use_rep = "X")

sc.tl.umap(eRegulon_gene_AUC)

sc.pl.umap(eRegulon_gene_AUC, color = "scRNA_counts:second_most_contributing_topic")

from scenicplus.RSS import (regulon_specificity_scores, plot_rss)

rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scRNA_counts:second_most_contributing_topic",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC"]
)

sc.pl.umap(eRegulon_gene_AUC, color = list(set([x for xs in [rss.loc[ct].sort_values()[0:2].index for ct in rss.index] for x in xs ])), save='_eregulon_enrichment_scores_genes.pdf')
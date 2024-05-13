from pycisTopic.cistopic_class import CistopicObject
from scipy.sparse import csr_matrix
from scipy.io import mmread
import gzip
from sklearn.preprocessing import binarize
import pandas as pd
import os
from pycisTopic.lda_models import run_cgs_models_mallet
from numpy.random import choice
from pycisTopic.lda_models import evaluate_models
import pickle
import numpy as np


"""
To build a CistopicObject:

class CistopicObject:
    
    cisTopic data class.

    :class:`CistopicObject` contains the cell by fragment matrices (stored as counts :attr:`fragment_matrix` and as binary accessibility :attr:`binary_matrix`),
    cell metadata :attr:`cell_data`, region metadata :attr:`region_data` and path/s to the fragments file/s :attr:`path_to_fragments`.

    LDA models from :class:`CisTopicLDAModel` can be stored :attr:`selected_model` as well as cell/region projections :attr:`projections` as a dictionary.

    Attributes
    ----------
    fragment_matrix: sparse.csr_matrix
        A matrix containing cell names as column names, regions as row names and fragment counts as values.
    binary_matrix: sparse.csr_matrix
        A matrix containing cell names as column names, regions as row names and whether regions as accessible (0: Not accessible; 1: Accessible) as values.
    cell_names: list
        A list containing cell names.
    region_names: list
        A list containing region names.
    cell_data: pd.DataFrame
        A data frame containing cell information, with cells as indexes and attributes as columns.
    region_data: pd.DataFrame
        A data frame containing region information, with region as indexes and attributes as columns.
    path_to_fragments: str or dict
        A list containing the paths to the fragments files used to generate the :class:`CistopicObject`.
    project: str
        Name of the cisTopic project.
"""

# Load and read matrix
print("Loading sparse matrix..")
fragment_matrix_gzip = gzip.open('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/matrices/monocyte/matrix.mtx.gz', 'rb')
coo_fragment_matrix = mmread(fragment_matrix_gzip)
fragment_matrix = csr_matrix(coo_fragment_matrix)
fragment_matrix_gzip.close()
print("Loaded sparse matrix!")

# Sample 5000 cells random
# Sample range based on columns (cells)
set.seed(639245)
sample_range = range(0,125985)
sample = choice(sample_range, 50000)

# Create dummy matrix
fragment_matrix_dummy = fragment_matrix[:,sample]

# Binarization sparse matrix 
print("Binarization of sparse matrix..")
fragment_matrix_binarized = binarize(fragment_matrix, threshold=0)
fragment_matrix_dummy_binarized = binarize(fragment_matrix_dummy, threshold=0)
print("Binarization done!")

print("Loading files for cisTopic object..")
# Read region names
region_names_file = gzip.open('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/matrices/monocyte/features.tsv.gz', 'rb')
region_names = pd.read_csv(region_names_file, sep='\t', header=None).iloc[:,0].to_list()
# region_names_dummy = [region_names[i] for i in sample]

# Read cell names
cell_names_file = gzip.open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/matrices/monocyte/barcodes.tsv.gz", 'rb')
cell_names = pd.read_csv(cell_names_file, sep='\t', header=None).iloc[:,0].to_list()
cell_names_dummy = [cell_names[i] for i in sample]

# Create dictionary with cell_names and index to later subset the fragments matrix as well
cell_names_index_dict = dict(zip(cell_names_dummy, range(0,len(cell_names_dummy))))
# Read region data
region_data = pd.read_csv("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/signac_peaks/output/mo_peaks_lane1to64_monocyte.bed", sep="\t")
region_data.index = region_names
# region_data_dummy = region_data[region_data.index.isin(region_names_dummy)]

# Fragment count cPeaks path
fragment_path = "/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/matrices/monocyte/"

# Cell metadata
cell_metadata = pd.read_csv("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_cellevel_metadata.tsv", sep="\t", dtype={'scrublet_doublet': str, 'filtered_condition': str})
# Sort barcodes
cell_metadata = cell_metadata.sort_values(by='barcode_lane')
# Subset cell metadata 
cell_metadata = cell_metadata[cell_metadata.barcode_lane.isin(cell_names)]
# Filter cell names
cell_names = cell_metadata.barcode_lane.to_list()
# Namees as rownames 
cell_metadata.index = cell_metadata.barcode_lane
cell_metadata_dummy = cell_metadata[cell_metadata.index.isin(cell_names_dummy)]
# Final filter cell names
cell_names_dummy = [i for i in cell_names_index_dict.keys() if i in cell_metadata_dummy.index]
# 4,459 cells left
# Check indices sampled cells that are left
sampled_remain = [cell_names_index_dict.get(i) for i in cell_names_dummy]
# Final filter fragment matrices
fragment_matrix_dummy_filtered = fragment_matrix_dummy[:,np.array(sampled_remain)]
fragment_matrix_dummy_binarized_filtered = fragment_matrix_dummy_binarized[:,np.array(sampled_remain)]

# Initialize cisTopic object
cistopic_obj = CistopicObject(
    fragment_matrix=fragment_matrix_dummy_filtered,
    binary_matrix=fragment_matrix_dummy_binarized_filtered,
    cell_names=cell_names_dummy,
    region_names=region_names,
    region_data=region_data,
    cell_data=cell_metadata_dummy,
    path_to_fragments=fragment_path,
    project="Monocytes_subset_4459_20_topics"
)
print("Finished cisTopic object!")
print(cistopic_obj)

# Run mallet for LDA
os.environ['MALLET_MEMORY'] = '200G'

# Mallet path
mallet_path = "/groups/umcg-franke-scrna/tmp03/users/umcg-mwvanderwerff/scenicplus/mallet/Mallet-202108/bin/mallet"

print("Runnign Mallet LDA..")
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[20],
    n_cpu=4,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path=os.environ["TMPDIR"],
    save_path="/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_50000_cells_20_topics",
    mallet_path=mallet_path,
)

print("Finished Mallet LDA")

model = evaluate_models(
    models,
    select_model = 20,
    return_model = True,
    save='figures/model_evaluation_mono_20_topics.pdf'
)
# Add model to cistopic object 
cistopic_obj.add_LDA_model(model)

# Write cistopic object with model to a pickle file
with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_50000_cells_20_topics/cistopic_obj.pkl", 'wb') as f:
   pickle.dump(cistopic_obj, f)

# Open cistopic object with model again is you start a new session
with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_50000_cells_20_topics/cistopic_obj.pkl", "rb") as f:
    cistopic_obj = pickle.load(f)

# Clustering and visualization
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)

# Identify clusters
find_clusters(
    cistopic_obj,
    target  = 'cell',
    k = 10,
    res = [0.6, 1.2, 3],
    prefix = 'pycisTopic_',
    scale = True,
    split_pattern = '_'
)

# Run UMAP
run_umap(
    cistopic_obj,
    target  = 'cell', scale=True)

# Plot topics for monocytes
plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['cell_type_lowerres_imputed', 'pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'],
    target='cell', num_columns=4,
    text_size=10,
    dot_size=5,
    save='figures/umap_20_topics_monocytes.pdf')

# Annotate with scRNA-seq data from Seurat
annot_dict = {}
for resolution in [0.6, 1.2, 3]:
    annot_dict[f"pycisTopic_leiden_10_{resolution}"] = {}
    for cluster in set(cistopic_obj.cell_data[f"pycisTopic_leiden_10_{resolution}"]):
        counts = cistopic_obj.cell_data.loc[
            cistopic_obj.cell_data.loc[cistopic_obj.cell_data[f"pycisTopic_leiden_10_{resolution}"] == cluster].index,
            "cell_type_lowerres_imputed"].value_counts()
        annot_dict[f"pycisTopic_leiden_10_{resolution}"][cluster] = f"{counts.index[counts.argmax()]}({cluster})"

for resolution in [0.6, 1.2, 3]:
    cistopic_obj.cell_data[f'pycisTopic_leiden_10_{resolution}'] = [
        annot_dict[f'pycisTopic_leiden_10_{resolution}'][x] for x in cistopic_obj.cell_data[f'pycisTopic_leiden_10_{resolution}'].tolist()
    ]

plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['cell_type_lowerres_imputed', 'pycisTopic_leiden_10_0.6', 'pycisTopic_leiden_10_1.2', 'pycisTopic_leiden_10_3'],
    target='cell', num_columns=2,
    text_size=10,
    dot_size=5,
    save='figures/umap_20_topics_monocytes_rna_annotated.pdf'
    )

# Continues plot
#plot_metadata(
#    cistopic_obj,
#    reduction_name='UMAP',
#    variables=['log10_unique_fragments_count', 'tss_enrichment', 'Doublet_scores_fragments', 'fraction_of_fragments_in_peaks'],
#    target='cell', num_columns=4,
#    text_size=10,
#    dot_size=5)

# Plot monocytes and color for topic contribution
plot_topic(
    cistopic_obj,
    reduction_name = 'UMAP',
    target = 'cell',
    num_columns=4,
    save='figures/umap_20_topics_monocytes_colored_by_topics.pdf'
)

# Heatmap of topics, does not work yet
cell_topic_heatmap(
    cistopic_obj,
    variables = ['cell_type_lowerres_imputed'],
    legend_loc_x = 1.0,
    legend_loc_y = -1.2,
    legend_dist_y = -1,
    figsize = (10, 10)
)

from pycisTopic.topic_binarization import binarize_topics

# Region bin topics plot binarization
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=4, save='figures/region_bin_plot_topics.pdf'
)

binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100,
    save='figures/binarized_cell_topics.pdf'
    )

"""
Following, we can compute the topic quality control metrics. These include:

    Number of assignments

    Topic coherence (Mimno et al., 2011): Measures to which extent high scoring regions in the topic are actually co-accessible in the original data. If it is low it indicates that the topic is rather random. The higher, the better is a topic.

    The marginal topic distribution: Indicates how much each topic contributes to the model. The higher, the better is a topic.

    The gini index: Value between 0 and 1, that indicates the specificity of topics (0: General, 1:Specific)

    If topics have been binarized, the number of regions/cells per topic will be added.

"""

from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img

topic_qc_metrics = compute_topic_metrics(cistopic_obj)

# Create dictionary with all QC figures
fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True, save='figures/coherence_vs_assignments_20_topics_5000_monocytes.pdf')
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True, save='figures/assignments_vs_cell_in_bin_20_topics_5000_monocytes.pdf')
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True, save='figures/coherence_vs_cell_in_bin_20_topics_5000_monocytes.pdf')
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True, save='figures/coherence_vs_regions_20_topics_5000_monocytes.pdf')
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True, save='figures/coherence_vs_marginal_dist_20_topics_5000_monocytes.pdf')
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True, save='figures/coherence_vs_gini_index_20_topics_5000_monocytes.pdf')

cistopic_obj.cell_data['cell_type_lowerres_imputed_mono'] = 'monocyte'

topic_annot = topic_annotation(
    cistopic_obj,
    annot_var='cell_type_lowerres_imputed',
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)

# DAR identification
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)

# Impute regions
imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)

# Normalize object
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True,
    save='figures/variable_regions_20_topics_5000_monocytes.pdf'
)

# Get the amount of variable regions
len(variable_regions)

# Find variable regions for stimulated and unstimulated monocytes
markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='filtered_condition_imputed',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=4,
    _temp_dir=os.environ["TMPDIR"],
    split_pattern = '_'
)

from pycisTopic.clust_vis import plot_imputed_features

plot_imputed_features(
    cistopic_obj,
    reduction_name='UMAP',
    imputed_data=imputed_acc_obj,
    features=[markers_dict[x].index.tolist()[0] for x in ['24hCa', 'UT']],
    scale=False,
    num_columns=4,
    save='figures/DARS_24hCa_UT_20_topics_5000_monocytes.pdf'
)

# Write cistopic object with model to a pickle file
with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_5000_cell_20_topics/cistopic_obj_newest.pkl", 'wb') as f:
   pickle.dump(cistopic_obj, f)

with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_5000_cell_20_topics/cistopic_obj_newest.pkl", "rb") as f:
    cistopic_obj = pickle.load(f)
# Plot stimulation status and other metrics
# Plot topics for monocytes

azimuth_annotation = pd.read_csv('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/cell_type_assignment/azimuth/10x_multiome_PBMCs/mo_azimuth_ct_10xmultiome.tsv', sep='\t')
azimuth_annotation = azimuth_annotation[azimuth_annotation['barcode'].isin(cistopic_obj.cell_data['barcode_lane'])]

azimuth_annotation = azimuth_annotation.sort_values(by=['barcode'])
cistopic_obj.cell_data = cistopic_obj.cell_data.sort_values(by=['barcode_lane'])
cistopic_obj.cell_data['predicted.mo_10x_cell_type'] = azimuth_annotation['predicted.mo_10x_cell_type'].to_list()

plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['filtered_condition_imputed','predicted.mo_10x_cell_type'],
    target='cell', num_columns=2,
    text_size=10,
    dot_size=5,
    save='figures/umap_condition_cell_sub_type.pdf')

binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100)

topic_annot = topic_annotation(
    cistopic_obj,
    annot_var='predicted.mo_10x_cell_type',
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)

# Add topic to which a cell contributes the most to the cell_data for each cell
cell_topics = cistopic_obj.selected_model.cell_topic.T
# Ignoring highest topic 16 for every cell, so getting second largest topic contribution
second_highest_topic = cell_topics.apply(lambda row: row.nlargest(2).index[-1],axis=1)
# Add to cell_data \
cistopic_obj.cell_data['second_most_contributing_topic'] = second_highest_topic

# Save updated cistopic object
with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_5000_cell_20_topics/cistopic_obj_newest.pkl", 'wb') as f:
   pickle.dump(cistopic_obj, f)

# Get list of all topics (second most contributing)
topics_dars_compare_with_topic_11 = cistopic_obj.cell_data.second_most_contributing_topic.unique().tolist()
# Remove topic 11 which we want to compare 
topics_dars_compare_with_topic_11.pop(7)
# Imput query DARs
input_query = [[['Topic11'], ['Topic7']], [['Topic11'], ['Topic10']], [['Topic11'], ['Topic8']],
               [['Topic11'], ['Topic6']], [['Topic11'], ['Topic19']], [['Topic11'], ['Topic20']],
               [['Topic11'], ['Topic17']], [['Topic11'], ['Topic15']], [['Topic11'], ['Topic18']],
               [['Topic11'], ['Topic3']], [['Topic11'], ['Topic12']], [['Topic11'], ['Topic5']],
               [['Topic11'], ['Topic13']]]

# Run DAR analysis
markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='second_most_contributing_topic',
    contrasts=input_query,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=4,
    _temp_dir=os.environ["TMPDIR"],
    split_pattern = '_'
)

from pycisTopic.clust_vis import plot_imputed_features

# Get non-empty results (some dicts were emtpy --> no DARs)
marker_dicts_not_empty = {k:v for (k,v) in markers_dict.items() if not v.empty}

# Get all keys of our markers_dict
comparisons = list(marker_dicts_not_empty.keys())

plot_imputed_features(
    cistopic_obj,
    reduction_name='UMAP',
    imputed_data=imputed_acc_obj,
    features=[markers_dict[x].index.tolist()[0] for x in comparisons],
    scale=False,
    num_columns=3,
    save='figures/DARs_topic_11_vs_others.pdf'
)

# Number of DARs
print("Number of DARs found:")
print("---------------------")
for x in marker_dicts_not_empty:
    print(f"  {x}: {len(marker_dicts_not_empty[x])}")

    
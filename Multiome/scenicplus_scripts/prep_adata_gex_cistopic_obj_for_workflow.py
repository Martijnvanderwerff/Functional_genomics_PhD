import scanpy as sc
import pandas as pd
import pickle

# Load adata from 10X format
adata = sc.read_10x_mtx(
    "/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/seurat_10x_format"
)

# Get cell data from metadata folder and subset for cells
cell_data = pd.read_csv("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/metadata/mo_cellevel_metadata.tsv", sep="\t")
cell_data_subset_keeps = cell_data['barcode_lane'].isin(adata.obs_names)
cell_data_subset = cell_data[cell_data_subset_keeps]

# Reorder subset dataset to match order of adata
df_adata = pd.DataFrame({'Cells' : adata.obs_names.to_list()})
cells_data_subset_reordered = pd.merge(df_adata, cell_data_subset,left_on='Cells',right_on='barcode_lane',how='outer')
cells_data_subset_reordered.index = cells_data_subset_reordered.Cells
adata.obs = cells_data_subset_reordered
# Needed for scenic+, crashed otherwise when pre-processing data
adata.raw = adata
# remove corrupted columns
del adata.obs['scrublet_doublet']
adata.write("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/adata_gex.h5ad")

###################################################################
# FURTHER PRE-PROCESS THE CISTOPIC OBJECT TO MATCH SCENIC+ WORKFLOW
###################################################################

# Remove index name cistopic_object
with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_5000_cell_20_topics/cistopic_obj.pkl", "rb") as f:
    cistopic_obj = pickle.load(f)

# Remove index name
cistopic_obj.cell_data.index.name=None

# Change region_names of cistopic object to chr:start-stop, otherwise error
cistopic_obj.region_names = [i.replace('-', ':', 1) for i in cistopic_obj.region_names]
cistopic_obj.region_data.index = [i.replace('-', ':', 1) for i in cistopic_obj.region_data.index]
cistopic_obj.region_data['name'] = [i.replace('-', ':', 1) for i in cistopic_obj.region_data['name']]
cistopic_obj.selected_model.topic_region.index = [i.replace('-', ':', 1) for i in cistopic_obj.selected_model.topic_region.index]

with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_5000_cell_20_topics/cistopic_obj.pkl", 'wb') as f:
   pickle.dump(cistopic_obj, f)
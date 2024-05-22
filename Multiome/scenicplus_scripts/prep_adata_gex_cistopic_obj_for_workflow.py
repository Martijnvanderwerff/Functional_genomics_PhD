import scanpy as sc
import pandas as pd
import pickle

# Load adata from 10X format
adata = sc.read_10x_mtx(
    "/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/mono_10000/seurat_10x_format"
)
# Get cell_data from our cistopic_obj
with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_10000_cells_20_topics/cistopic_obj.pkl", "rb") as f:
    cistopic_obj = pickle.load(f)
cell_data = cistopic_obj.cell_data

# Add observations dataframe to adata.obs
adata.obs = pd.merge(adata.obs, cell_data, right_index=True, left_index=True)
# Needed for scenic+, crashed otherwise when pre-processing data
adata.raw = adata

adata.write("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/mono_10000/adata_gex.h5ad")

# CODE BELOW IS NOT NEEDED ANYMORE

###################################################################
# FURTHER PRE-PROCESS THE CISTOPIC OBJECT TO MATCH SCENIC+ WORKFLOW
###################################################################

# Remove index name cistopic_object
# with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_5000_cell_20_topics/cistopic_obj.pkl", "rb") as f:
#     cistopic_obj = pickle.load(f)

# # Remove index name
# cistopic_obj.cell_data.index.name=None

# Change region_names of cistopic object to chr:start-stop, otherwise error
#cistopic_obj.region_names = [i.replace('-', ':', 1) for i in cistopic_obj.region_names]
#cistopic_obj.region_data.index = [i.replace('-', ':', 1) for i in cistopic_obj.region_data.index]
#cistopic_obj.region_data['name'] = [i.replace('-', ':', 1) for i in cistopic_obj.region_data['name']]
#cistopic_obj.selected_model.topic_region.index = [i.replace('-', ':', 1) for i in cistopic_obj.selected_model.topic_region.index]
#cistopic_obj.selected_model.region_topic.index = [i.replace('-', ':', 1) for i in cistopic_obj.selected_model.region_topic.index]


# with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_5000_cell_20_topics/cistopic_obj.pkl", 'wb') as f:
#    pickle.dump(cistopic_obj, f)
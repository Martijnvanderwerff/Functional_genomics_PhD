import scanpy as sc
import pandas as pd
import pickle

# Load adata from 10X format
adata = sc.read_10x_mtx(
    "/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/mono_all/seurat_10x_format"
)
# Get cell_data from our cistopic_obj
# with open("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_10000_cells_20_topics/cistopic_obj.pkl", "rb") as f:
#     cistopic_obj = pickle.load(f)
# cell_data = cistopic_obj.cell_data

# Add observations dataframe to adata.obs
# adata.obs = pd.merge(adata.obs, cell_data, right_index=True, left_index=True)
# Needed for scenic+, crashed otherwise when pre-processing data
adata.raw = adata

adata.write("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/mono_all/adata_gex.h5ad")

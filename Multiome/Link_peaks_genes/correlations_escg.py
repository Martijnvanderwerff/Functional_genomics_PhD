
from scipy.sparse import csr_matrix
from scipy.io import mmread
from scipy.io import mmwrite
import gzip
import pandas as pd
import numpy as np

fragment_matrix_gzip = gzip.open('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/matrices/monocyte/matrix.mtx.gz', 'rb')
coo_fragment_matrix = mmread(fragment_matrix_gzip)
fragment_matrix = csr_matrix(coo_fragment_matrix)
fragment_matrix_gzip.close()

region_names_file = gzip.open('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/matrices/monocyte/features.tsv.gz', 'rb')
region_names = pd.read_csv(region_names_file, sep='\t', header=None).iloc[:,0].to_list()
# Get index and position in region_names list
chrom_6_regions = [(i,j) for i,j in enumerate(region_names) if j.startswith('chr6')]
chrom_6_regions_indices = [i[0] for i in chrom_6_regions]
chrom_6_matrix = fragment_matrix[chrom_6_regions_indices,:]

# Write chrom6 matrix to file
mmwrite('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/chr6_regions.mtx', chrom_6_matrix)

# Get matrix peak regions
chrom_6_matrix = mmread('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/chr6_regions.mtx')
# Keep regions if 0.001 cells have a region
# Get number of non-zero values per columns
non_zero_per_row = chrom_6_matrix.getnnz(axis=1)
# Calculate fraction compared to number of cells
non_zero_per_row_fraction = non_zero_per_row / chrom_6_matrix.shape[1]
# Filter for cells >= 0.001
regions_pass_filter = np.argwhere(non_zero_per_row_fraction >= 0.001).ravel()
# Filter chrom 6 matrix on regions where 0.001 cells have that region
chrom_6_matrix_filtered = chrom_6_matrix[regions_pass_filter,:]
# Write filtered matrix 
mmwrite('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/chr6_regions_filtered.mtx', chrom_6_matrix_filtered)
# Get final region names 
chrom_6_regions_filtered = [chrom_6_regions[i] for i in regions_pass_filter]
# Write filterd region set to txt file
with open('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/chr6_regions_filtered.txt', 'w') as f:
    for line in chrom_6_regions_filtered:
        f.write(f"{line[1]}\n")

# Read cell names from monocyte matrix
cell_names_atac = pd.read_csv("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/matrices/monocyte/barcodes.tsv.gz", sep="\t", header=None).iloc[:,0].to_list()
# Read RNA count matrix 
rna_matrix = mmread('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/seurat_10x_format/matrix.mtx')
rna_matrix = csr_matrix(rna_matrix)
rna_cells = pd.read_csv("/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/seurat_object/seurat_10x_format/barcodes.tsv", sep="\t", header=None).iloc[:,0].to_list()

# Get matching cell indices for the RNA data, in order to subset
rna_cell_in_atac_cells = []
for i, x in enumerate(rna_cells):
    if x in cell_names_atac:
        rna_cell_in_atac_cells.append(i)
        print(i)

# Filter RNA matrix so that ATAC and RNA cells match
rna_matrix_filtered = rna_matrix[:, rna_cell_in_atac_cells]
# Write RNA matrix with matching cells to file
mmwrite('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/rna_genes/monocytes_filtered.mtx', rna_matrix_filtered)
# Write cells to file of all monoyctes 
with open('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/region_gene_linking/monocytes/rna_genes/monocyte_cells_filtered.txt', 'w') as f:
    for line in cell_names_atac:
        f.write(f"{line}\n")
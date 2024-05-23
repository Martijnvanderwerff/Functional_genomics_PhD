import mudata
import pandas as pd

# Load data
scplus_mdata = mudata.read('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/scenic_workflow/test_10000_monocytes_20_topics/outs/scplusmdata.h5mu')
directed_e_regulons = scplus_mdata.uns["direct_e_regulon_metadata"]
extended_e_regulons = scplus_mdata.uns["extended_e_regulon_metadata"]

directed_e_regulons.to_csv('direct_e_regulons.csv', index=None)
extended_e_regulons.to_csv('extended_e_regulons.csv', index=None)

# Load Topic5 DARs 
dars_topic5 = pd.read_csv('/groups/umcg-franke-scrna/tmp03/projects/multiome/ongoing/scenic_plus/pycistopic/monocytes_10000_cells_20_topics/region_sets/DARS_topics_vs_others/Topic5.bed', sep='\t')
dars_topic5['Region'] = 
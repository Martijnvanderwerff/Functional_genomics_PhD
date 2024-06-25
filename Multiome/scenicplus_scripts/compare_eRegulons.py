import pandas as pd

# Open direct eRegulons for 10 topics and 10 topics models monocytes 

mono_10_topics = pd.read_csv('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scenicplus_pipeline/monocytes/model_10_topics/Snakemake/outs/eRegulon_direct.tsv', sep='\t')
mono_20_topics = pd.read_csv('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scenicplus_pipeline/monocytes/model_20_topics/Snakemake/outs/eRegulon_direct.tsv', sep='\t')

mono_10_topics[['Region', 'Gene']].isin(mono_20_topics[['Region', 'Gene']]).value_counts()

s1 = pd.merge(mono_10_topics, mono_20_topics, how='inner', on=['Region', 'Gene', 'TF', 'eRegulon_name'])

# Identify DARs that are in CRE's (percentage)

# 10 topics all mono

# Import Module 
import os 
import numpy as np
# Folder Path 
path_model_10 = "/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/pycistopic/monocytes_all_cells/region_sets/model_10_topics/DARS_between_topics"
path_model_20 = "/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/pycistopic/monocytes_all_cells/region_sets/model_20_topics/DARS_topics_vs_others"

# Read text File 
def read_text_file(file_path): 
  
  with open(file_path, 'r') as f: 
      lines = [line.rstrip() for line in f]
      lines = [sub.replace('\t', '-') for sub in lines]
      lines = [sub.replace('-', ':', 1) for sub in lines]
  return lines

# Add regions to list 
def comparse_dars(model_name, path, cres):
  dars_pycistopic_all= []
  # iterate through all file 
  for file in os.listdir(path=path): 
    # Check whether file is in text format or not 
    if file.endswith(".bed"): 
      file_path = f"{path}/{file}"
      # call read text file function 
      dars_pycistopic_all.append(read_text_file(file_path))
  dars_pycistopic_all = np.concatenate(dars_pycistopic_all).tolist()
  regions_cre = cres['Region'].tolist()
  dars_in_cre = set(regions_cre).intersection(set(dars_pycistopic_all))
  print(f'Percent DARs in CREs {model_name} topics: {len(dars_in_cre) / len(dars_pycistopic_all) * 100}')

comparse_dars('model 10', path_model_10, mono_10_topics)
comparse_dars('model 20', path_model_20, mono_20_topics)
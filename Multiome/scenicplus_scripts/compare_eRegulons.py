import pandas as pd
import seaborn as sns
import colorcet as cc
import os 
import numpy as np
from venn import venn
import matplotlib.pyplot as plt
# Open direct eRegulons for 10 topics and 10 topics models monocytes 

mono_10_topics = pd.read_csv('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scenicplus_pipeline/monocytes/model_10_topics/Snakemake/outs/eRegulon_direct.tsv', sep='\t')
mono_20_topics = pd.read_csv('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scenicplus_pipeline/monocytes/model_20_topics/Snakemake/outs/eRegulon_direct.tsv', sep='\t')
mono_10k_20_topics = pd.read_csv('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scenicplus_pipeline/monocytes/model_10k_mono_20_topics/eRegulon_direct.tsv', sep='\t')

mono_10_topics[['Region', 'Gene']].isin(mono_20_topics[['Region', 'Gene']]).value_counts()

s1 = pd.merge(mono_10_topics, mono_20_topics, how='inner', on=['Region', 'Gene', 'TF', 'eRegulon_name'])

# Identify DARs that are in CRE's (percentage)

# 10 topics all mono

# Import Module 

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

# Check overlap DARs 10k, all mono 10/20 topics
mono_10k_dars = read_text_file('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/pycistopic/monocytes_10000_cells_20_topics/region_sets/DARS_topics_vs_others/all_dars_10k_mono.bed')
mono_all_10_topics_dars = read_text_file('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/pycistopic/monocytes_all_cells/region_sets/model_10_topics/DARS_between_topics/all_dars.bed')
mono_all_20_topics_dars = read_text_file('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/pycistopic/monocytes_all_cells/region_sets/model_20_topics/DARS_topics_vs_others/all_dars.bed')

datas = {
    "Mono 10k": set(mono_10k_dars),
    "Mono all 10 topics": set(mono_all_10_topics_dars),
    "Mono all 20 topics": set(mono_all_20_topics_dars)
}

venn(datas).get_figure().savefig('overlap_dars.png')

# Fetch common DARs for all three methods
common_dars = set(mono_10k_dars).intersection(mono_all_10_topics_dars, mono_all_20_topics_dars)

def common_dars_check(model_name, cres, dars):
  regions_cre = cres['Region'].tolist()
  dars_in_cre = set(regions_cre).intersection(set(dars))
  print(f'Common DARs found in CREs for {model_name}: {len(dars_in_cre)}/{len(common_dars)}')

common_dars_check('all monocytes 10 topics', mono_10_topics, common_dars)
common_dars_check('all monocytes 20 topics', mono_20_topics, common_dars)
common_dars_check('10k monocytes 20 topics', mono_10k_20_topics, common_dars)

# Get unique eRegulon_directs for each dataframe model 
unique_mono_10_topics = mono_10_topics[(~mono_10_topics['eRegulon_name'].isin(mono_20_topics['eRegulon_name'])) & (~mono_10_topics['eRegulon_name'].isin(mono_10k_20_topics['eRegulon_name']))]
unique_mono_20_topics = mono_20_topics[(~mono_20_topics['eRegulon_name'].isin(mono_10_topics['eRegulon_name'])) & (~mono_20_topics['eRegulon_name'].isin(mono_10k_20_topics['eRegulon_name']))]
unique_mono_10k_20_topics = mono_10k_20_topics[(~mono_10k_20_topics['eRegulon_name'].isin(mono_10_topics['eRegulon_name'])) & (~mono_10k_20_topics['eRegulon_name'].isin(mono_20_topics['eRegulon_name']))]

# Save DARs that are in CREs as bed files, so we have subset files 
def common_dars_check(cres, dars):
  regions_cre = cres['Region'].tolist()
  dars_in_cre = set(regions_cre).intersection(set(dars))
  return dars_in_cre

# Change : and - to \t again
mono_10k_dars_in_eRegulons = [sub.replace('-', '\t').replace(':', '\t') for sub in common_dars_check(mono_10k_20_topics, mono_10k_dars)]
mono_all_10_topics_dars_in_eRegulons = [sub.replace('-', '\t').replace(':', '\t') for sub in common_dars_check(mono_10_topics, mono_all_10_topics_dars)]
mono_all_20_topics_dars_in_eRegulons = [sub.replace('-', '\t').replace(':', '\t') for sub in common_dars_check(mono_20_topics, mono_all_20_topics_dars)]

with open('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scripts/regions_in_eRegulons/mono_10k.bed', 'w') as f:
    for line in mono_10k_dars_in_eRegulons:
        f.write(f"{line}\n")
with open('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scripts/regions_in_eRegulons/mono_10_topics.bed', 'w') as f:
    for line in mono_all_10_topics_dars_in_eRegulons:
        f.write(f"{line}\n")
with open('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scripts/regions_in_eRegulons/mono_20_topics.bed', 'w') as f:
    for line in mono_all_20_topics_dars_in_eRegulons:
        f.write(f"{line}\n")


# Plot bed files with annotation as upsetplot for DARs that are in eRegulons

def plot_annotations_bars(annotation_file, output_file):
  data = pd.read_csv(annotation_file, sep='\t', header=None)
  data = data[3].value_counts()
  data.index = [item.split("/")[-1] for item in data.index]
  data.index = [item.split("=")[-1] for item in data.index]
  data.index = [item.split(".")[0] for item in data.index]
  data = data.sort_index(ascending=True)
  # Plot annotations
  fig, ax = plt.subplots()
  bars = ax.barh(data.index.values, data.values)
  ax.bar_label(bars)
  plt.savefig(output_file, bbox_inches='tight')
  plt.close()
  return data 


data_10k = plot_annotations_bars(annotation_file='/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scripts/eRegulon_region_annotations/mono_10k_20_topics.bed', output_file='/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/figures/barplot_annotations_mono_10k.pdf')
data_all_10_topics = plot_annotations_bars(annotation_file='/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scripts/eRegulon_region_annotations/mono_all_10_topics.bed', output_file='/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/figures/barplot_annotations_all_mono_10_topics.pdf')
data_all_20_topics = plot_annotations_bars(annotation_file='/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scripts/eRegulon_region_annotations/mono_all_20_topics.bed', output_file='/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/figures/barplot_annotations_all_mono_20_topics.pdf')

fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(10, 11))
bars1 = ax1.barh(data_10k.index.values, data_10k.values)
ax1.bar_label(bars1)
bars2 = ax2.barh(data_all_10_topics.index.values, data_all_10_topics.values)
ax2.bar_label(bars2)
bars3 = ax3.barh(data_all_20_topics.index.values, data_all_20_topics.values)
ax3.bar_label(bars3)
ax1.title.set_text(f'10k monocytes subset, {data_10k.sum()} regions')
ax2.title.set_text(f'All monocytes 10 topics, {data_all_10_topics.sum()} regions')
ax3.title.set_text(f'All monocytes 20 topics, {data_all_20_topics.sum()} regions')
plt.suptitle('Annotated regions used in eRegulons')
plt.savefig('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/figures/barplot_annotations_grouped.pdf', bbox_inches='tight')
plt.close()


def plot_annotations_pie(annotation_file, output_file):
  data = pd.read_csv(annotation_file, sep='\t', header=None)
  data = data[3].value_counts()
  data.index = [item.split("/")[-1] for item in data.index]
  data.index = [item.split("=")[-1] for item in data.index]
  data.index = [item.split(".")[0] for item in data.index]
  data = data.sort_index(ascending=True)
  fig = plt.figure()
  patches, texts = plt.pie(data)
  percent = 100.*data/data.sum()
  labels = labels = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(data.index, percent)]
  sort_legend = True
  if sort_legend:
      patches, labels, dummy =  zip(*sorted(zip(patches, labels, data),
                                          key=lambda x: x[2],
                                          reverse=True))
  plt.legend(patches, labels, bbox_to_anchor=(-0.1, 1.), fontsize=7, loc='center')  
  plt.savefig(output_file, bbox_inches='tight')
  plt.close()
  return data

data_10k = plot_annotations_bars(annotation_file='/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scripts/eRegulon_region_annotations/mono_10k_20_topics.bed', output_file='/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/figures/piechart_annotations_mono_10k.pdf')
data_all_10_topics = plot_annotations_bars(annotation_file='/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scripts/eRegulon_region_annotations/mono_all_10_topics.bed', output_file='/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/figures/piechart_annotations_all_mono_10_topics.pdf')
data_all_20_topics = plot_annotations_bars(annotation_file='/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/scripts/eRegulon_region_annotations/mono_all_20_topics.bed', output_file='/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/figures/piechart_annotations_all_mono_20_topics.pdf')

# Get housekeeping data
data_10k_housekeeping = data_10k[['NOThousekeeping', 'IShousekeeping']]
data_all_10_topics_housekeeping = data_all_10_topics[['NOThousekeeping', 'IShousekeeping']]
data_all_20_topics_housekeeping = data_all_20_topics[['NOThousekeeping', 'IShousekeeping']]

# Get not housekeeping data
data_10k_not_housekeeping = data_10k.drop(labels=['NA','NOThousekeeping', 'IShousekeeping', 'NAhousekeeping'])
data_all_10_topics_not_housekeeping = data_all_10_topics.drop(labels=['NA', 'NOThousekeeping', 'IShousekeeping', 'NAhousekeeping'])
data_all_20_topics_not_housekeeping = data_all_20_topics.drop(labels=['NA', 'NOThousekeeping', 'IShousekeeping', 'NAhousekeeping'])

# Plot housekeeping regions
fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3)
colors = sns.color_palette('deep')
ax1.pie(data_10k_housekeeping,
        autopct="%1.2f%%",
        labels=data_10k_housekeeping.index,
        startangle=90,
        explode=tuple([0.1] * len(data_10k_housekeeping.index)),
        colors=colors,
        textprops={'fontsize': 9})

ax2.pie(data_all_10_topics_housekeeping,
        autopct="%1.2f%%",
        labels=data_all_10_topics_housekeeping.index,
        startangle=90,
        explode=tuple([0.1] * len(data_all_10_topics_housekeeping.index)),
        colors=colors,
        textprops={'fontsize': 9})

ax3.pie(data_all_20_topics_housekeeping,
        autopct="%1.2f%%",
        labels=data_all_20_topics_housekeeping.index,
        startangle=90,
        explode=tuple([0.1] * len(data_all_20_topics_housekeeping.index)),
        colors=colors,
        textprops={'fontsize': 9})

ax1.set_title(f'10k monocytes subset, {data_10k_housekeeping.sum()} regions\n\nNA={data_10k["NA"]}, NAhousekeeping={data_10k["NAhousekeeping"]}', fontsize=9)
ax2.set_title(f'All monocytes 10 topics, {data_all_10_topics_housekeeping.sum()} regions\n\nNA={data_all_10_topics["NA"]}, NAhousekeeping={data_all_10_topics["NAhousekeeping"]}', fontsize=9)
ax3.set_title(f'All monocytes 20 topics, {data_all_20_topics_housekeeping.sum()} regions\n\nNA={data_all_20_topics["NA"]}, NAhousekeeping={data_all_20_topics["NAhousekeeping"]}', fontsize=9)
plt.suptitle('Housekeeping classification regions used in eRegulons', weight='bold', y=0.75)
plt.savefig('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/figures/pie_annotations_grouped_housekeeping.pdf', bbox_inches='tight')
plt.close()

# Plot non-housekeeping regions
color_dict_not_housekeeping = { label: color  for label, color in zip(data_10k_not_housekeeping.index, sns.color_palette(cc.glasbey, n_colors=16))}

fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3)
patches, texts = ax1.pie(data_10k_not_housekeeping,
        colors=[color_dict_not_housekeeping[key] for key in data_10k_not_housekeeping.index])
percent = 100.*data_10k_not_housekeeping/data_10k_not_housekeeping.sum()
labels = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(data_10k_not_housekeeping.index, percent)]
patches, labels, dummy =  zip(*sorted(zip(patches, labels, data_10k_not_housekeeping),
                                      key=lambda x: x[2],
                                      reverse=True))
ax1.legend(patches, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=1, fontsize=10)

patches_2, texts_2 = ax2.pie(data_all_10_topics_not_housekeeping,
        colors=[color_dict_not_housekeeping[key] for key in data_all_10_topics_not_housekeeping.index])
percent_2 = 100.*data_all_10_topics_not_housekeeping/data_all_10_topics_not_housekeeping.sum()
labels_2 = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(data_all_10_topics_not_housekeeping.index, percent_2)]
patches_2, labels_2, dummy_2 =  zip(*sorted(zip(patches_2, labels_2, data_all_10_topics_not_housekeeping),
                                      key=lambda x: x[2],
                                      reverse=True))
ax2.legend(patches_2, labels_2, loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=1, fontsize=10)


patches_3, texts_3 = ax3.pie(data_all_20_topics_not_housekeeping,
        colors=[color_dict_not_housekeeping[key] for key in data_all_20_topics_not_housekeeping.index])
percent_3 = 100.*data_all_20_topics_not_housekeeping/data_all_20_topics_not_housekeeping.sum()
labels_3 = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(data_all_20_topics_not_housekeeping.index, percent_3)]
patches_3, labels_3, dummy_3 =  zip(*sorted(zip(patches_3, labels_3, data_all_20_topics_not_housekeeping),
                                      key=lambda x: x[2],
                                      reverse=True))
ax3.legend(patches_3, labels_3, loc='upper center', bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True, ncol=1, fontsize=10)

ax1.set_title(f'10k monocytes subset, {data_10k_not_housekeeping.sum()} regions', fontsize=11)
ax2.set_title(f'All monocytes 10 topics, {data_all_10_topics_not_housekeeping.sum()} regions', fontsize=11)
ax3.set_title(f'All monocytes 20 topics, {data_all_20_topics_not_housekeeping.sum()} regions', fontsize=11)
plt.suptitle('Non-housekeeping classification regions used in eRegulons', weight='bold', y=0.75)
plt.savefig('/scratch/hb-functionalgenomics/projects/multiome/ongoing/scenicplus_workdir/figures/pie_annotations_grouped_not_housekeeping.pdf', bbox_inches='tight')
plt.close()

module load utility
source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/pycistopic_polars_3_10_env
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

import os
import pickle
import pandas as pd

#infer gene activity
import pyranges as pr
from pycisTopic.gene_activity import get_gene_activity

num_core = 1

cancer = "oligo"
base_dir = '/mnt/beegfs/home/alyne/stagein/scenicplus_'+cancer
os.chdir(base_dir)

out_dir = '/mnt/beegfs/home/alyne/stageout/pycisTopic_output_'+cancer
fig_dir = out_dir+"/figures_pycistopic_"+cancer
pct_obj_dir = out_dir+"/pct_objects"

#get size of chromosomes from UCSC
chromsizes = pd.read_table(
    "../required_files/chrom_sizes.txt",
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)
chromsizes.head()
chromsizes = pr.PyRanges(chromsizes[["Chromosome", "Start", "End"]])


#load pycistopic object
with open(os.path.join(pct_obj_dir, "cistopic_obj.pkl"), "rb") as f:
    cistopic_obj = pickle.load(f)



#binarise topic regions and cell topic distributions
#various different methods available
#first binarise topic region distributions
from pycisTopic.topic_binarization import binarize_topics

region_bin_topics_otsu = binarize_topics(
    cistopic_obj.selected_model.topic_region, 
    cistopic_obj.selected_model.topic_region.index,
    method='otsu'
)
#also do ntop, might be required for later part of pipeline
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj.selected_model.topic_region, 
    cistopic_obj.selected_model.topic_region.index,
    method='ntop', ntop = 3000
)

#then binarize cells
binarized_cell_topic = binarize_topics(
    cistopic_obj.selected_model.cell_topic.T,
    cistopic_obj.selected_model.cell_topic.T.index,
    method='li')


#annotate topics with cell types by computing proportion of cells in each group assigned the annotation
#and comparing to overall percentage

# Build dict in the expected format
binarized_dict = {
    f"Topic{i+1}": pd.Series(
        binarized_cell_topic[1][i],  # scores
        index=binarized_cell_topic[0][i]  # matching cell IDs
    )
    for i in range(len(binarized_cell_topic[0]))
}

# Now run annotation
from pycisTopic.topic_qc import topic_annotation

topic_annot = topic_annotation(
    cistopic_obj,
    annot_var="cell_type",
    binarized_cell_topic=binarized_dict
)


import numpy as np
from pycisTopic.cistopic_class import *
from pycisTopic.diff_features import *


topic_region_ = cistopic_obj.selected_model.topic_region.to_numpy(dtype=np.float32)
cell_topic_ = cistopic_obj.selected_model.cell_topic.to_numpy(dtype=np.float32)
regions = cistopic_obj.selected_model.topic_region.index.tolist()


(region_names_to_keep, per_region_means_on_normalized_imputed_acc, per_region_dispersions_on_normalized_imputed_acc,) = calculate_per_region_mean_and_dispersion_on_normalized_imputed_acc(
          region_topic=topic_region_,
          cell_topic=cell_topic_,
          region_names=regions,
          scale_factor1 = 10**6,
          scale_factor2 = 10**4,
          regions_chunk_size=20000,)


np.savez_compressed(out_dir+"/ct_output/impute_acc_per_region_mean.npz", arr=per_region_means_on_normalized_imputed_acc)
np.savez_compressed(out_dir+"/ct_output/impute_acc_per_region_dispersion.npz", arr=per_region_dispersions_on_normalized_imputed_acc)
np.savez_compressed(out_dir+"/ct_output/region_idx_to_keep.npz", arr=region_names_to_keep)
pd.DataFrame(region_names_to_keep, columns=["regions"]).to_csv(out_dir+"/ct_output/region_names_to_keep.csv", index=False)


## Find highly variable features.
var_features = find_highly_variable_features(
          features=region_names_to_keep,
          per_region_means_on_normalized_imputed_acc=per_region_means_on_normalized_imputed_acc,
          per_region_dispersions_on_normalized_imputed_acc=per_region_dispersions_on_normalized_imputed_acc,
          min_disp = 0.025,
          min_mean = 0.01,
          max_disp = float("inf"),
          max_mean = 3,
          n_bins = 20,
          n_top_features = None,
          plot = fig_dir + "/variable_features.png")

pd.DataFrame(var_features, columns=["regions"]).to_csv(out_dir+"/ct_output/var_features.csv", index=False)


## use old impute_accessibility function with variable features only
## note that this step may still require quite some resources
## this is only an intermediate solution and will be updated soon 
imputed_acc_obj = impute_accessibility(
          cistopic_obj,
          selected_cells=None,
          selected_regions=var_features,
          scale_factor=10**6
)

with open(out_dir+"/ct_output/imputed_cistopic_obj.pkl", 'wb') as f:
          pickle.dump(imputed_acc_obj, f)
# #load imputed pycistopic object
# with open(out_dir+"/ct_output/imputed_cistopic_obj.pkl", "rb") as f:
#     imputed_acc_obj = pickle.load(f)

## calculate DARs
markers_dict = find_diff_features(
          cistopic_obj,
          imputed_acc_obj, 
          variable='cell_type', var_features=var_features, contrasts = None,
          adjpval_thr=0.05, log2fc_thr=np.log2(1.5),
          n_cpu=1, _temp_dir='tmp_dir'
)

with open(out_dir+"/ct_output/DARs.pkl", 'wb') as f:
          pickle.dump(markers_dict, f)
# #load markers dict
# with open(out_dir+"/ct_output/DARs.pkl", "rb") as f:
#     markers_dict = pickle.load(f)

from pycisTopic.clust_vis import plot_imputed_features
#plot region accessibility on the cell-topic UMAP. 
#For example, let’s check how the best DARs for some cell types look.
features = sum(
    [markers_dict[x].index.tolist()[:2] for x in
     ['Macrophage','GABAergic','Glutamatergic','Astro','Oligo',
      'OPC-like','intermediate','astro-like']],
    []
)

plot_imputed_features(
    cistopic_obj,
    reduction_name='UMAP',
    imputed_data=imputed_acc_obj,
    features=features,
    scale=False,
    num_columns=5,
    save=fig_dir + "/umaps_best_DARs_per_cell_type.png"
)

#plot top 9 for each cell type in separate figure (NB v slow)
cell_types = ['Macrophage','GABAergic','Glutamatergic','Astro','Oligo',
    'OPC-like','intermediate','astro-like']
for ct in cell_types:
    plot_imputed_features(cistopic_obj, 'UMAP', imputed_acc_obj,
                           markers_dict[ct].index.tolist()[:9],
                           scale=False, num_columns=3,
                           save=f"{fig_dir}/UMAP_top9_DARs_{ct}.png")


#plot aggregated score for top 100 DARs for each cell type
agg = {}
for ct in cell_types:
    # get top 100 DARs for this cell type
    feats = markers_dict[ct].index.tolist()[:100]
    # subset the imputed object for these features
    sub = imputed_acc_obj.subset(features=feats, copy=True)
    # get the matrix (mtx or mat depending on your version)
    mat = sub.mtx
    # compute mean accessibility across these features for each cell
    agg[ct] = np.asarray(mat).mean(axis=0)

# convert the dictionary into a DataFrame
agg_df = pd.DataFrame(agg, index=imputed_acc_obj.cell_names)
# attach aggregated signals to cistopic_obj metadata
for ct in cell_types:
    cistopic_obj.cell_data[f"{ct}_top100_DAR_mean"] = agg_df[ct].values

from pycisTopic.clust_vis import plot_metadata
# plot all aggregated signals on UMAP
plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=[f"{ct}_top100_DAR_mean" for ct in cell_types],
    num_columns=4,
    save=f"{fig_dir}/UMAP_aggregated_top100_DARs_per_celltype.png"
)



#save region sets
os.makedirs(os.path.join(out_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "DARs_cell_type"), exist_ok = True)

from pycisTopic.utils import region_names_to_coordinates

for topic in range(len(region_bin_topics_otsu[0])):
    region_names_to_coordinates(
        region_bin_topics_otsu[0][topic]
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for topic in range(len(region_bin_topics_top_3k[0])):
    region_names_to_coordinates(
        region_bin_topics_top_3k[0][topic]
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for cell_type in markers_dict:
    region_names_to_coordinates(
        markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )


pickle_path = os.path.join(out_dir, "region_sets", "Topics_otsu", "region_bin_topics_otsu.pkl")
with open(pickle_path, "wb") as f:
    pickle.dump(region_bin_topics_otsu, f)

pickle_path = os.path.join(out_dir, "region_sets", "Topics_otsu", "region_bin_topics_top_3k.pkl")
with open(pickle_path, "wb") as f:
    pickle.dump(region_bin_topics_top_3k, f)

pickle_path = os.path.join(out_dir, "region_sets", "Topics_otsu", "markers_dict.pkl")
with open(pickle_path, "wb") as f:
    pickle.dump(markers_dict, f)



#infer gene activity
from pycisTopic.gene_activity import get_gene_activity

pr_annotation = pd.read_table(
        os.path.join("..", "required_files", "tss.bed")
    ).rename(
        {"Name": "Gene", "# Chromosome": "Chromosome"}, axis = 1)
pr_annotation["Transcription_Start_Site"] = pr_annotation["Start"]
pr_annotation = pr.PyRanges(pr_annotation)


from pycisTopic.diff_features import (impute_accessibility,find_highly_variable_features,find_diff_features)
import numpy as np


#impute full accessibility
imputed_acc_obj = impute_accessibility(
          cistopic_obj,
          selected_cells=None,
          selected_regions=None,
          scale_factor=10**6
)

#infer gene accessibility
gene_act, weigths = get_gene_activity(
    imputed_acc_obj,
    pr_annotation,
    chromsizes,
    use_gene_boundaries=True, # Whether to use the whole search space or stop when encountering another gene
    upstream=[1000, 100000], # Search space upstream. The minimum means that even if there is a gene right next to it
                             # these bp will be taken (1kbp here)
    downstream=[1000,100000], # Search space downstream
    distance_weight=True, # Whether to add a distance weight (an exponential function, the weight will decrease with distance)
    decay_rate=1, # Exponent for the distance exponential funciton (the higher the faster will be the decrease)
    extend_gene_body_upstream=10000, # Number of bp upstream immune to the distance weight (their value will be maximum for
                          #this weight)
    extend_gene_body_downstream=500, # Number of bp downstream immune to the distance weight
    gene_size_weight=False, # Whether to add a weights based on the length of the gene
    gene_size_scale_factor='median', # Dividend to calculate the gene size weigth. Default is the median value of all genes
                          #in the genome
    remove_promoters=False, # Whether to remove promoters when computing gene activity scores
    average_scores=True, # Whether to divide by the total number of region assigned to a gene when calculating the gene
                          #activity score
    scale_factor=1, # Value to multiply for the final gene activity matrix
    extend_tss=[10,10], # Space to consider a promoter
    gini_weight = True, # Whether to add a gini index weigth. The more unique the region is, the higher this weight will be
    return_weights= True, # Whether to return the final weights
    project='Gene_activity') # Project name for the gene activity object

pickle.dump(
    gene_act,
    open(os.path.join(out_dir, "gene_act_obj.pkl"), "wb")
)


# #also do differentially accessible genes
# DAG_markers_dict = find_diff_features(
#     cistopic_obj, gene_act, variable='cell_type', var_features=None,
#     contrasts=None, adjpval_thr=0.05,
#     log2fc_thr=np.log2(1.5), n_cpu=num_core, split_pattern = '-')
# 

from pycisTopic.clust_vis import plot_imputed_features

plot_imputed_features(
    cistopic_obj, reduction_name='UMAP', imputed_data=gene_act,
    features=['JUN', 'FOS', 'DIAPH3', "ADGRV1", 'KCNN3', 'TNC', 'MYT1L', 'ROBO1', 'ROBO2','FGF12', 'LHFPL3', 'TNR'], 
    scale=True, num_columns=4, save = fig_dir + f"/umaps_imputed_features_{cancer}_markers.pdf"
)


















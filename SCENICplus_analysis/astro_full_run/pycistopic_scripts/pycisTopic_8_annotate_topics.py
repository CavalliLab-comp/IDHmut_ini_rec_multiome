module load utility
source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/pycistopic_polars_3_10_env

import os
import pickle
import pandas as pd

cancer = "astro"
base_dir = "/mnt/beegfs/home/alyne"
os.chdir(base_dir)


#directories to store output
out_dir = base_dir+'/stageout/pycisTopic_output_'+cancer
fig_dir = out_dir+"/figures_pycistopic_"+cancer
pct_obj_dir = out_dir+"/pct_objects"

#load in models
topic_counts = [30,40,50,60,70,80]

# Container for models
models = []
for k in topic_counts:
    model_path = out_dir+f"/models/tm.{k}_topics.model.pkl"
    with open(model_path, "rb") as f:
        model = pickle.load(f)
    models.append(model)


#evaluate models using different methods
from pycisTopic.lda_models import evaluate_models
model = evaluate_models(
    models, select_model = 60, return_model = True, plot = False,
    save = fig_dir+"/evaluate_topic_models_"+cancer+".pdf"
)


with open(os.path.join(pct_obj_dir, "cistopic_obj.pkl"), "rb") as f:
    cistopic_obj = pickle.load(f)
    

cistopic_obj.add_LDA_model(model)


from pycisTopic.clust_vis import (
    find_clusters, run_umap, run_tsne, plot_metadata, plot_topic, cell_topic_heatmap
)

ress = [0.3,0.5,0.8]
#run clustering and umap using topic contribution
find_clusters(
    cistopic_obj, target = 'cell', k = 10, res = ress,
    prefix = 'pycisTopic_', scale = True, split_pattern = '-'
)
run_umap(cistopic_obj, target = 'cell', scale=True)

import seaborn as sns
categories = cistopic_obj.cell_data['cell_type'].astype(str).str.strip().unique()
palette = sns.color_palette("hls", n_colors=len(categories))

color_dict = {
    'cell_type': dict(zip(categories, palette))
}
plot_metadata(
    cistopic_obj,reduction_name='UMAP', show_label = False, show_legend = True,
    variables=['cell_type','sample','pycisTopic_leiden_10_0.3','pycisTopic_leiden_10_0.5','pycisTopic_leiden_10_0.8'], 
    target='cell', num_columns=3, text_size=10, dot_size=1, 
    save = fig_dir+"/umaps_"+cancer+"_cell_type.pdf", color_dictionary = color_dict)



#annotate clusters based on scRNA-seq annotations
annot_dict = {}
for resolution in ress:
    annot_dict[f"pycisTopic_leiden_10_{resolution}"] = {}
    for cluster in set(cistopic_obj.cell_data[f"pycisTopic_leiden_10_{resolution}"]):
        counts = cistopic_obj.cell_data.loc[
            cistopic_obj.cell_data.loc[cistopic_obj.cell_data[f"pycisTopic_leiden_10_{resolution}"] == cluster].index,
            "cell_type"].value_counts()
        annot_dict[f"pycisTopic_leiden_10_{resolution}"][cluster] = f"{counts.index[counts.argmax()]}({cluster})"


for resolution in ress:
    cistopic_obj.cell_data[f'pycisTopic_leiden_10_{resolution}'] = [
        annot_dict[f'pycisTopic_leiden_10_{resolution}'][x] for x in cistopic_obj.cell_data[f'pycisTopic_leiden_10_{resolution}'].tolist()
    ]


#plot cluster annotation
plot_metadata(
    cistopic_obj, reduction_name='UMAP', variables=['cell_type', 'pycisTopic_leiden_10_0.3','pycisTopic_leiden_10_0.5','pycisTopic_leiden_10_0.8'],
    target='cell', num_columns=2, text_size=7, dot_size=1, color_dictionary = color_dict,
    save = fig_dir+"/umaps_"+cancer+"_cluster_annotation.pdf")


#plot QC metrics
plot_metadata(
    cistopic_obj, reduction_name='UMAP',
    variables=['log10_unique_fragments_count', 'tss_enrichment', 'fraction_of_fragments_in_peaks'],
    target='cell', num_columns=3, text_size=10, dot_size=1,
    save = fig_dir+"/umaps_"+cancer+"_QC_metrics.pdf")


#umap of topic contribution
plot_topic(
    cistopic_obj, reduction_name = 'UMAP', target = 'cell', num_columns=8,
    save = fig_dir+"/umaps_"+cancer+"_topic_scores.pdf"
)


#heatmap of topic contribtion
import matplotlib.pyplot as plt
figsize = (20, 15)
cell_topic_heatmap(
    cistopic_obj, variables = ['cell_type'], scale = False, legend_loc_x = 0.4,
    legend_loc_y = -1.2, legend_dist_y = -1, figsize = figsize,
    color_dictionary = color_dict, save = None
)
# Adjust y-axis ticks and labels
ax = plt.gca()  # Get current axis
ax.set_yticks(range(len(your_labels)))  # Adjust this based on your data
ax.set_yticklabels(your_labels)  # Replace with your actual labels

# Additional adjustments
plt.subplots_adjust(left=0.25)
plt.tick_params(axis='y', labelsize=4)  # Adjust font size
plt.yticks(rotation=45)  # Rotate labels for better visibility

plt.savefig(fig_dir + f"/heatmap_{cancer}_topic_scores.pdf", bbox_inches="tight")
plt.close()


#save again
pickle.dump(
    cistopic_obj,
    open(os.path.join(pct_obj_dir, "cistopic_obj.pkl"), "wb")
)


# compute topic quality control metrics
from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img

topic_qc_metrics = compute_topic_metrics(cistopic_obj)

fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
# fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
# fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
# fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)

fig=plt.figure(figsize=(12, 10))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(1, 3, i)
    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1

plt.subplots_adjust(wspace=0, hspace=-0.5)
fig.savefig(fig_dir + f"/topic_{cancer}_QC_metrics.pdf", bbox_inches="tight", dpi = 300)



#get matrix of topic contributions so can identify which topics go with which cell type
import pandas as pd

cell_topic_df = pd.DataFrame(
    cistopic_obj.selected_model.cell_topic
)
#add cell types/states
cell_topic_df = cell_topic_df.T
cell_topic_df['cell_type'] = cistopic_obj.cell_data['cell_type']

#compute mean topic score per cell type
celltype_topic_scores = (cell_topic_df
    .groupby('cell_type')
    .mean()
)
#find top 3 topics
top3_topics = {ct: (celltype_topic_scores.loc[ct]
        .sort_values(ascending=False)
        .head(3)
        .index
        .tolist()
    )
    for ct in celltype_topic_scores.index
}
top3_topics_df = pd.DataFrame(top3_topics)
top3_topics_df.to_csv(fig_dir+"/top3_topics_per_celltype_"+cancer+".tsv",header=True,sep="\t")



#export matrix for use in ArchR to analyse ATAC only
from scipy import sparse
# get region names
regions = cistopic_obj.region_names
#convert to bed format
bed = pd.DataFrame([r.replace(":", "-").split("-") for r in regions],columns=["chr","start","end"])
bed.to_csv(out_dir+"/pycistopic_peaks_"+cancer+".bed", sep="\t", header=False, index=False)
# # matrix
# mat = cistopic_obj.fragment_matrix
# # save sparse matrix
# sparse.save_npz(out_dir+"/cistopic_peak_matrix_"+cancer+".npz", mat)
# # save peaks
# pd.Series(cistopic_obj.region_names).to_csv(out_dir+"/peaks_"+cancer+".txt", index=False, header=False)
# # save cells
# pd.Series(cistopic_obj.cell_names).to_csv(out_dir+"/cells_"+cancer+".txt", index=False, header=False)
# # also save cell names in ArchR format
# cells = cistopic_obj.cell_names
# barcodes = [c.split("-astro")[0] for c in cells]
# samples = [c.split("___")[1] for c in cells]
# archr_cells = [f"{s}#{b}" for s,b in zip(samples,barcodes)]
# pd.Series(archr_cells).to_csv(out_dir+"/archr_cells_"+cancer+".txt", index=False, header=False)
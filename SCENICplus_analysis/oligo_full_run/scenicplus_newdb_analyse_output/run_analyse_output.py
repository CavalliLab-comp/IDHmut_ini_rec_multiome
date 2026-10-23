module load utility
source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
# conda create --prefix /mnt/beegfs/home/alyne/persistent/scenicplus_after_env -c conda-forge python=3.11.8
conda activate /mnt/beegfs/home/alyne/persistent/scenicplus_after_env
##pip install .
##pip install scanpy --no-deps
##pip install anndata numpy scipy pandas matplotlib seaborn h5py
####pip install "pandas==2.0.3"

import os
cancer = "oligo"

out_dir = "/mnt/beegfs/home/alyne/stageout/scenicplus_op_newdb_filt_"+cancer
fig_dir = out_dir+"/figures_scenicplus_newdb_"+cancer
os.makedirs(fig_dir, exist_ok = True)
os.chdir(fig_dir)

import mudata
scplus_mdata = mudata.read(out_dir+"/outs/scplusmdata.h5mu")
scplus_mdata.uns["direct_e_regulon_metadata"]
scplus_mdata.uns["extended_e_regulon_metadata"]


# umap based on regulon AUC scores
#eregulon dimensionality reduction
import scanpy as sc
import anndata
eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,
)

eRegulon_gene_AUC.obs = scplus_mdata.obs.loc[eRegulon_gene_AUC.obs_names]
sc.pp.neighbors(eRegulon_gene_AUC, use_rep = "X")
sc.tl.umap(eRegulon_gene_AUC)
sc.pl.umap(eRegulon_gene_AUC, color = "scATAC_counts:cell_type",
  save = "_TFactivity.png")


#find regulons with highest specificity per cell type
#eregulon specificity score
from scenicplus.RSS import (regulon_specificity_scores, plot_rss)

rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scATAC_counts:cell_type",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC"]
)

plot_rss(
    data_matrix = rss, top_n = 4, num_columns = 4, save = "figures/regulon_specificity.png"
)
rss.to_csv("figures/rss_"+cancer+".csv", sep = ';')


#plot eregulon enrichment score for top n for each cell type
#number to plot for each cell type
n_top = 2
#find top regulons and make titles
colors = []
titles = []
for ct in rss.index:
    top_regs = (
        rss.loc[ct]
           .sort_values(ascending=False)
           .head(n_top)
           .index
           .tolist()
    )
    for reg in top_regs:
        colors.append(reg)
        titles.append(f"{ct} – {reg}")

sc.pl.umap(
    eRegulon_gene_AUC, color=colors, title=titles, save="_enriched_regulons.png"
)



#heatmap dotplot

from scenicplus.plotting.dotplot import heatmap_dotplot

heatmap_dotplot(
    scplus_mudata = scplus_mdata,
    color_modality = "direct_gene_based_AUC",
    size_modality = "direct_region_based_AUC",
    group_variable = "scATAC_counts:cell_type",
    eRegulon_metadata_key = "direct_e_regulon_metadata",
    color_feature_key = "Gene_signature_name",
    size_feature_key = "Region_signature_name",
    feature_name_key = "eRegulon_name",
    sort_data_by = "direct_gene_based_AUC",
    orientation = "horizontal",
    figsize = (25, 20), save = "figures/heatmap_dotplot.pdf"
)


#write out specificity scores and activity scores for 
import pandas as pd
import numpy as np
# # Make sure it's dense
# X = eRegulon_gene_AUC.X.toarray() if not isinstance(eRegulon_gene_AUC.X, np.ndarray) else eRegulon_gene_AUC.X


#write gene_based AUC scores to file
df = pd.DataFrame(scplus_mdata["direct_gene_based_AUC"].X, 
	index = scplus_mdata["direct_gene_based_AUC"].obs_names, 
	columns = scplus_mdata["direct_gene_based_AUC"].var_names)
df.to_csv("direct_gene_based_AUC_"+cancer+".tsv.gz", index=True, sep="\t", compression="gzip")

df_ex = pd.DataFrame(scplus_mdata["extended_gene_based_AUC"].X, 
	index = scplus_mdata["extended_gene_based_AUC"].obs_names, 
	columns = scplus_mdata["extended_gene_based_AUC"].var_names)
df_ex.to_csv("extended_gene_based_AUC_"+cancer+".tsv.gz", index=True, sep="\t", compression="gzip")

df = pd.DataFrame(scplus_mdata["direct_region_based_AUC"].X, 
	index = scplus_mdata["direct_region_based_AUC"].obs_names, 
	columns = scplus_mdata["direct_region_based_AUC"].var_names)
df.to_csv("direct_region_based_AUC_"+cancer+".tsv.gz", index=True, sep="\t", compression="gzip")

df_ex = pd.DataFrame(scplus_mdata["extended_region_based_AUC"].X, 
	index = scplus_mdata["extended_region_based_AUC"].obs_names, 
	columns = scplus_mdata["extended_region_based_AUC"].var_names)
df_ex.to_csv("extended_region_based_AUC_"+cancer+".tsv.gz", index=True, sep="\t", compression="gzip")



df = scplus_mdata.uns["direct_e_regulon_metadata"]
df.to_csv("direct_regulon_md_"+cancer+".tsv.gz", index=True, sep="\t", compression="gzip")

df_ex = scplus_mdata.uns["extended_e_regulon_metadata"]
df_ex.to_csv("extended_regulon_md_"+cancer+".tsv.gz", index=True, sep="\t", compression="gzip")







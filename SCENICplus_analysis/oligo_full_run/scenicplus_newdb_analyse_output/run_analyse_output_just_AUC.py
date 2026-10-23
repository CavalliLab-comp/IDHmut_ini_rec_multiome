
import os
import mudata
import pandas as pd


cancer = "astro"
fig_dir = "/Users/alyne/Library/CloudStorage/Dropbox/Curie/Curie_StCloud/anne-marie/10x_RNA_ATAC_oligo_analysis/SCENICplus_soroban/astro_full_run/figures_scenicplus_astro"

aucell_dir_mdata = mudata.read("/Users/alyne/Library/CloudStorage/Dropbox/Curie/Curie_StCloud/anne-marie/"+
"10x_RNA_ATAC_oligo_analysis/SCENICplus_soroban/"+cancer+"_full_run/outs_op_newdb_final/AUCell_direct.h5mu")
aucell_ext_mdata = mudata.read("/Users/alyne/Library/CloudStorage/Dropbox/Curie/Curie_StCloud/anne-marie/"+
"10x_RNA_ATAC_oligo_analysis/SCENICplus_soroban/"+cancer+"_full_run/outs_op_newdb_final/AUCell_extended.h5mu")



import scanpy as sc
import anndata
gene_AUC = anndata.concat(
    [aucell_dir_mdata["Gene_based"], aucell_ext_mdata["Gene_based"]],
    axis = 1,
)
#adapt names to match anotation fle
gene_AUC.obs_names = gene_AUC.obs_names.str.split('___').str[0]


#write output
import scipy.sparse as sp

# Convert X to dense if it's sparse
if sp.issparse(gene_AUC.X):
    X = gene_AUC.X.toarray()
else:
    X = gene_AUC.X

# Make DataFrame
df = pd.DataFrame(X,index=gene_AUC.obs_names,columns=gene_AUC.var_names)

# Save to TSV
df.to_csv(fig_dir+"/../outs_op_newdb_final/gene_AUC_sp_"+cancer+".tsv.gz", sep="\t")



cell_annot = pd.read_csv(fig_dir+'/../outs_op_newdb_final/cell_annotations.tsv.gz', compression='gzip', sep='\t',index_col=0)

gene_AUC.obs = gene_AUC.obs.join(cell_annot, how='left')


os.chdir(fig_dir)

sc.pp.neighbors(gene_AUC, use_rep = "X")
sc.tl.umap(gene_AUC)
sc.pl.umap(gene_AUC, color = "cell_type", save = "_TFactivity.png")





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

# Cells as rows, genes as columns
df = pd.DataFrame(eRegulon_gene_AUC.X, index=eRegulon_gene_AUC.obs_names, 
  columns=eRegulon_gene_AUC.var_names)
df.to_csv("figures/eRegulon_gene_AUC_astro_test.csv", sep = ';')

rss.to_csv("figures/rss_astro_test.csv", sep = ';')


#write gene_based AUC scores to file
df = pd.DataFrame(scplus_mdata["direct_gene_based_AUC"].X, 
	index = scplus_mdata["direct_gene_based_AUC"].obs_names, 
	columns = scplus_mdata["direct_gene_based_AUC"].var_names)
df.to_csv("direct_gene_based_AUC_"+cancer+".tsv.gz", index=True, sep="\t", compression="gzip")

df_ex = pd.DataFrame(scplus_mdata["extended_gene_based_AUC"].X, 
	index = scplus_mdata["extended_gene_based_AUC"].obs_names, 
	columns = scplus_mdata["extended_gene_based_AUC"].var_names)
df_ex.to_csv("extended_gene_based_AUC_"+cancer+".tsv.gz", index=True, sep="\t", compression="gzip")




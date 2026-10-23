
module load utility
source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/scanpy_env

import anndata as ad
import os
import pandas as pd
from scipy import sparse


cancer = "oligo"
wd = "/mnt/beegfs/home/alyne/stagein/scenicplus_"+cancer+"/RNA"

os.chdir(wd)

raw_counts = pd.read_csv("raw_counts_scenicplus_"+cancer+".tsv.gz", sep="\t", index_col=0)
norm_expr = pd.read_csv("SCT_data_scenicplus_"+cancer+".tsv.gz", sep="\t", index_col=0)

raw_counts_sparse = sparse.csr_matrix(raw_counts.values)
norm_expr_sparse = sparse.csr_matrix(norm_expr.values)

adata = ad.AnnData(X=norm_expr_sparse)   # .X expects cells × genes, has already been transposed
adata.var_names = norm_expr.columns
adata.obs_names = norm_expr.index

adata.raw = ad.AnnData(X=raw_counts_sparse, var=pd.DataFrame(index=raw_counts.columns))

adata.write("adata_object_scenicplus_"+cancer+".h5ad")





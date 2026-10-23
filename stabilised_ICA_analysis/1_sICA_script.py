import scanpy as sc
import os
os.environ["NUMBA_DISABLE_CACHING"] = "1"
import pandas as pd
# from sica.singlecell import ica
from sicaomics.singlecell import ica
# import skdim
import numpy as np
from random import sample
import sys

samp = sys.argv[1]

wd = "/data/users/alyne/gliomult/svn/analyse/results/analysis_"+samp+"/sICA_"+samp
# wd = "/data/users/alyne/gliomult/svn/analyse/results/analysis_"+samp+"/sICA_filt_"+samp

os.chdir(wd)

adata = sc.read("/data/tmp/seurat_obj_aml/oligo_malig_"+samp+".tsv.gz",
delimiter='\t',first_column_names = True)
# adata = sc.read("/data/tmp/seurat_obj_filt/oligo_malig_"+samp+".tsv.gz",
# delimiter='\t',first_column_names = True)

adata.X -= adata.X.mean(axis=0)
sICA = ica(adata, observations = 'genes', n_components = 25, n_runs = 100,
		  return_info = True, copy = True, algorithm = 'picard_ext')

Metagenes = pd.DataFrame(sICA.varm['sica_metagenes'])
Metagenes.index = adata.var.index
Metagenes.to_csv('metagenes_'+samp+'.tsv.gz', index = True,
				  header = True, sep = '\t', compression='gzip')

Metasamples = pd.DataFrame(sICA.obsm['sica_metasamples'])
Metasamples.index = adata.obs.index
Metasamples.to_csv('metasamples_'+samp+'.tsv.gz', index = True,
				  header = True, sep = '\t', compression='gzip')

stability = pd.DataFrame(sICA.uns['sica']['stability_indexes'])
stability.to_csv('stability_'+samp+'.tsv.gz', index = True,
				  header = True, sep = '\t', compression='gzip')


#'fastica_par', 'fastica_def', 'picard_fastica', 'picard', 'picard_ext', 'picard_orth'

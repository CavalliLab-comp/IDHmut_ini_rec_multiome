


module load utility
source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/pycistopic_polars_3_10_env

#needs quite a lot of ram to make ct_objects (~400GB)

#number of cores available
num_core = 1

import os
import pickle
import pandas as pd

cancer = "oligo"
base_dir = '/mnt/beegfs/home/alyne/stagein/scenicplus_'+cancer
os.chdir(base_dir)

#read in list of sample names
#REMEMBER USE ONLY SAMPLES WITH RNA AND ATAC
samps = pd.read_csv("sample_list.txt", header=None)   # no header in file
sample_names = samps[0].tolist()

#make directories to store output
out_dir = '/mnt/beegfs/home/alyne/stageout/pycisTopic_output_'+cancer

fig_dir = out_dir+"/figures_pycistopic_"+cancer

# Create the fragment dictionary
fragments_dict = {
    sample: os.path.join(base_dir, f"analysis_{sample}", "CR_output", "atac_fragments.tsv.gz")
    for sample in sample_names
}

#read in cell annotation
import pandas as pd

cell_data = pd.read_table("cell_annotations.tsv.gz", index_col = 0)
cell_data.head()



#read in tsv file of cells passisng both QC filters
df = pd.read_csv("RNA/cells_passing_both_QC_"+cancer+".tsv.gz", sep="\t", compression="gzip")
sample_id_to_barcodes_passing_filters = {col: df[col].dropna().tolist() for col in df.columns}



#create cistopic object
#count matrix of fragments in consensus peaks for each cell

path_to_regions = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed")
path_to_blacklist = "../required_files/hg38-blacklist.v2.bed.gz"
pycistopic_qc_output_dir = out_dir+"/qc"
pct_obj_dir = out_dir+"/pct_objects"
os.makedirs(pct_obj_dir, exist_ok = True)

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl



#Step 1 make cistopic obj and save for each sample
for sample_id in fragments_dict:
    sample_metrics = pl.read_parquet(
        os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
    ).to_pandas().set_index("CB").loc[ sample_id_to_barcodes_passing_filters[sample_id] ]
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,
        path_to_blacklist = path_to_blacklist,
        metrics = sample_metrics,
        valid_bc = sample_id_to_barcodes_passing_filters[sample_id],
        n_cpu = num_core,
        project = sample_id,
        split_pattern = '-'
    )
    # Save to disk to avoid keeping all in memory
    with open(os.path.join(pct_obj_dir, f"{sample_id}_cistopic_obj.pkl"), "wb") as f:
        pickle.dump(cistopic_obj, f)



# for sample_id in list(fragments_dict.keys())[11:]:
#     sample_metrics = pl.read_parquet(
#         os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
#     ).to_pandas().set_index("CB").loc[ sample_id_to_barcodes_passing_filters[sample_id] ]
#     cistopic_obj = create_cistopic_object_from_fragments(
#         path_to_fragments = fragments_dict[sample_id],
#         path_to_regions = path_to_regions,
#         path_to_blacklist = path_to_blacklist,
#         metrics = sample_metrics,
#         valid_bc = sample_id_to_barcodes_passing_filters[sample_id],
#         n_cpu = num_core,
#         project = sample_id,
#         split_pattern = '-'
#     )
#     # Save to disk to avoid keeping all in memory
#     with open(os.path.join(pct_obj_dir, f"{sample_id}_cistopic_obj.pkl"), "wb") as f:
#         pickle.dump(cistopic_obj, f)
# 






import pickle
# Step 2: Load back and merge when ready
loaded_objs = []
for sample_id in fragments_dict:
    with open(os.path.join(pct_obj_dir, f"{sample_id}_cistopic_obj.pkl"), "rb") as f:
        loaded_objs.append(pickle.load(f))


#add cell metadata to each object
for obj in loaded_objs:
	obj.add_cell_data(cell_data, split_pattern='-')


#there seem to be some cells without annotation so remove them
for obj in loaded_objs:
	print(obj.cell_data['cell_type'].isna().sum())
	# Keep only cells with a non-missing cell_type
	cells_to_keep = obj.cell_data['cell_type'].notna()
	filtered_cells = obj.cell_data.index[cells_to_keep]
	# Subset the CistopicObject
	obj = obj.subset(filtered_cells)



#merge list to get one object
from pycisTopic.cistopic_class import merge
cistopic_obj = merge(loaded_objs)


pickle.dump(
    cistopic_obj,
    open(os.path.join(pct_obj_dir, "cistopic_obj.pkl"), "wb")
)


#write accessibility matrix to file
from scipy import sparse
from scipy.io import mmwrite
# get binary matrix
binary_matrix = cistopic_obj.binary_matrix
# write to file
binary_matrix_path=out_dir+"/ct_output"
os.makedirs(binary_matrix_path, exist_ok = True)
mmwrite(os.path.join(binary_matrix_path,"binary_accessibility_matrix.mtx"),binary_matrix)


# get region names and cell barcodes
regions = cistopic_obj.region_names
cells = cistopic_obj.cell_names

# write to files

regions_path=binary_matrix_path+'/region_names.txt'
with open(regions_path, 'w') as f:
    for name in regions:
        f.write(f"{name}\n")

cells_path=binary_matrix_path+'/cell_barcodes.txt'
with open(cells_path, 'w') as f:
    for name in cells:
        f.write(f"{name}\n")
        


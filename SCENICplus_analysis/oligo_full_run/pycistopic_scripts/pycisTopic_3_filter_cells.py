source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/pycistopic_polars_3_10_env


##############
### QC should just have been run in terminal
### for each sample
##############


import os
import pickle
import pandas as pd

cancer = "oligo"
base_dir = '/mnt/beegfs/home/alyne/stagein/scenicplus_'+cancer
os.chdir(base_dir)

#read in list of sample names
samps = pd.read_csv("sample_list.txt", header=None)   # no header in file
sample_names = samps[0].tolist()

#number of cores available
num_core = 10

#make directories to store output
out_dir = '/mnt/beegfs/home/alyne/stageout/pycisTopic_output_'+cancer
os.makedirs(out_dir, exist_ok = True)

fig_dir = out_dir+"/figures_pycistopic_"+cancer
os.makedirs(fig_dir, exist_ok = True)

# Create the fragment dictionary
fragments_dict = {
    sample: os.path.join(base_dir, f"analysis_{sample}", "CR_output", "atac_fragments.tsv.gz")
    for sample in sample_names
}

#read in cell annotation
import pandas as pd

# plot some qc figures
from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats
import matplotlib.pyplot as plt


for sample_id in fragments_dict:
    fig = plot_sample_stats( sample_id = sample_id, pycistopic_qc_output_dir = out_dir+"/qc",
    	save = fig_dir+"/QC_stats_"+sample_id+".pdf")


# fig = plot_sample_stats( sample_id = sample_id, pycistopic_qc_output_dir = out_dir+"/qc",
#     	save = fig_dir+"/QC_stats_"+sample_id+".pdf")

from pycisTopic.qc import get_barcodes_passing_qc_for_sample

sample_id_to_barcodes_passing_filters = {}
sample_id_to_thresholds = {}
for sample_id in fragments_dict:
    (
        sample_id_to_barcodes_passing_filters[sample_id],
        sample_id_to_thresholds[sample_id]
    ) = get_barcodes_passing_qc_for_sample(
            sample_id = sample_id,
            pycistopic_qc_output_dir = out_dir+"/qc",
            unique_fragments_threshold = None, # use automatic thresholding
            tss_enrichment_threshold = None, # use automatic thresholding
            frip_threshold = 0,
            use_automatic_thresholds = True,
    )


for sample_id in fragments_dict:
    fig = plot_barcode_stats(
        sample_id = sample_id,
        pycistopic_qc_output_dir = out_dir+"/qc",
        bc_passing_filters = sample_id_to_barcodes_passing_filters[sample_id],
        detailed_title = False,
        **sample_id_to_thresholds[sample_id],
        save = fig_dir+"/QC_barcode_stats_"+sample_id+".pdf"
    )


#write list of cells passing ATAC thresholds
import csv
fieldnames = sample_id_to_barcodes_passing_filters.keys()
csv_filename = out_dir+"/cells_passing_qc_ATAC.csv"

# Find the maximum list length
max_len = max(len(v) for v in sample_id_to_barcodes_passing_filters.values())
# Pad lists with "NaN"
padded_data = {k: v + ["NaN"] * (max_len - len(v)) for k, v in sample_id_to_barcodes_passing_filters.items()}

# Write to CSV
with open(csv_filename, "w", newline="") as f:
    writer = csv.writer(f, delimiter=";")
    # Write header
    writer.writerow(padded_data.keys())
    # Write rows
    for i in range(max_len):
        row = [padded_data[k][i] for k in padded_data]
        writer.writerow(row)



#########################
###### Now go to RNA, Plot overlap showing cell types
###### Make anndata object for input to scenicplus
#########################


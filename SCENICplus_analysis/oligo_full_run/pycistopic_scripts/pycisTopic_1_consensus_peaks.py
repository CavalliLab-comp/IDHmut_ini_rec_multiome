#script to be run once for all samples
#call consensus peaks for each cell type/state


module load utility
source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/pycistopic_polars_3_10_env

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

#number of cores available
#think this part maybe needs quite a few cores but can't get it to work
num_core = 10

#make directories to store output
out_dir = '/mnt/beegfs/home/alyne/stageout/pycisTopic_output_'+cancer
os.makedirs(out_dir, exist_ok = True)

fig_dir = out_dir+"/figures_scenicplus_"+cancer
os.makedirs(fig_dir, exist_ok = True)

# Create the fragment dictionary
fragments_dict = {
    sample: os.path.join(base_dir, f"analysis_{sample}", "CR_output", "atac_fragments.tsv.gz")
    for sample in sample_names
}

#read in cell annotation
import pandas as pd

cell_data = pd.read_table("cell_annotations.tsv.gz", index_col = 0)
cell_data.head()


#get size of chromosomes from UCSC
chromsizes = pd.read_table(
    "../required_files/chrom_sizes.txt",
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)
chromsizes.head()


#generate pseudobulk ATAC profiles
os.makedirs(os.path.join(out_dir, "consensus_peak_calling"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"), 
	exist_ok = True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"), 
	exist_ok = True)


from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
bw_paths, bed_paths = export_pseudobulk(
    input_data = cell_data, variable = "cell_type",
    sample_id_col = "sample", chromsizes = chromsizes,
    bed_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    bigwig_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"),
    path_to_fragments = fragments_dict, n_cpu = num_core,
    normalize_bigwig = True, temp_dir = "/tmp", split_pattern = "-"
)


#write paths of bed files to disk
with open(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv"), "wt") as f:
    for v in bw_paths:
        _ = f.write(f"{v}\t{bw_paths[v]}\n")


with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        _ = f.write(f"{v}\t{bed_paths[v]}\n")



#get paths to bed files if starting from here
import glob

tsv_folder = os.path.join(out_dir, "consensus_peak_calling", "pseudobulk_bed_files")
bed_paths = {}

for gz_path in glob.glob(os.path.join(tsv_folder, "*.tsv.gz")):
    sample = os.path.basename(gz_path).replace(".tsv.gz", "")
    bed_paths[sample] = gz_path

print("Pseudobulk files:", bed_paths)




#use MACS to call peaks for each pseudobulk
from pycisTopic.pseudobulk_peak_calling import peak_calling
macs_path = "macs2"

os.makedirs(os.path.join(out_dir, "consensus_peak_calling/MACS"), exist_ok = True)

narrow_peak_dict = peak_calling(
    macs_path = macs_path, bed_paths = bed_paths,
    outdir = os.path.join(out_dir, "consensus_peak_calling/MACS"),
    genome_size = 'hs', n_cpu = num_core, input_format = 'BEDPE',
    shift = 73, ext_size = 146, keep_dup = 'all', q_value = 0.05
)
 
        
#get set of consensus peaks
from pycisTopic.iterative_peak_calling import get_consensus_peaks
# Other param
peak_half_width=250
path_to_blacklist="../required_files/hg38-blacklist.v2.bed.gz"
# Get consensus peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes,
    path_to_blacklist = path_to_blacklist)


consensus_peaks.to_bed(
    path = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed"),
    keep = True, compression = 'infer', chain = False)


#QC for single cell ATAC-seq samples
#most of this to be run in terminal



# #make TSS anotation file in terminal
# have already run this so just transfer file to cluster
# awk 'BEGIN{OFS="\t"} 
#     $3=="transcript" {
#         split($0,a,";")
#         gene_name="NA"; gene_type="NA"
#         for(i in a){
#             if(a[i] ~ /gene_name/){split(a[i],b," "); gene_name=b[2]; gsub(/"/,"",gene_name)}
#             if(a[i] ~ /gene_type/){split(a[i],b," "); gene_type=b[2]; gsub(/"/,"",gene_type)}
#         }
#         if(gene_type=="protein_coding"){
#             tss=($7=="+")?$4-1:$5
#             print $1,tss,tss+1,gene_name,".",$7,gene_type
#         }
#     }' <(gzip -dc input_data/gencode.v49.annotation.gtf.gz) > input_data/tss.bed
# 


##############
#### Next run QC via terminal
#### pycisTopic_2_compute_qc.sh
##############




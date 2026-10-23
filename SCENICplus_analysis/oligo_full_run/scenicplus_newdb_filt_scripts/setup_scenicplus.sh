module load utility
source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/scenicplus_env
### 


cancer="oligo"

cd /mnt/beegfs/home/alyne/stageout
mkdir scenicplus_op_newdb_filt_"$cancer"

cd scenicplus_op_newdb_filt_"$cancer"

mkdir -p scplus_pipeline
mkdir outs
mkdir tmp
scenicplus init_snakemake --out_dir scplus_pipeline

### copy yaml file locally and complete input paths
### add full output path
### add chromsizes.tsv and gene_annotation.tsv to outs folder





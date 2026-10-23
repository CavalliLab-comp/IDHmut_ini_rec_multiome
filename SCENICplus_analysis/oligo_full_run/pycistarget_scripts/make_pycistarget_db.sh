#!/bin/bash
#SBATCH --job-name=ct_db
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=500gb
#SBATCH --time=48:00:00
#SBATCH --partition=dev
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err
#SBATCH --mail-type=END
#SBATCH --mail-user=anne-marie.lyne@curie.fr

source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/create_cistarget_db_env

cancer="oligo"
#### Variables
genome_fasta="/mnt/beegfs/common/annotations/pipelines_CDR/Human/hg38/genome/hg38.fa"
region_bed="/mnt/beegfs/home/alyne/stageout/pycisTopic_output_${cancer}/consensus_peak_calling/consensus_regions.bed"
region_fasta="/mnt/beegfs/home/alyne/stagein/required_files/consensus_regions_${cancer}.fa"
database_dir="/mnt/beegfs/home/alyne/stageout/pycistarget_db_${cancer}_cellstates_filt"
mkdir $database_dir
database_suffix="${database_dir}/db_${cancer}_cellstates_filt"
path_to_motif_collection="/mnt/beegfs/home/alyne/stagein/aertslab_motif_colleciton/v10nr_clust_public/singletons"

# ls $path_to_motif_collection > /mnt/beegfs/home/alyne/stagein/aertslab_motif_colleciton/v10nr_clust_public/motifs.txt
motif_list='/mnt/beegfs/home/alyne/stagein/aertslab_motif_colleciton/v10nr_clust_public/motifs.txt'
n_cpu='9'


#### Get fasta sequences
bedtools getfasta -fi ${genome_fasta} -bed ${region_bed} > ${region_fasta}


#### Set ${create_cistarget_databases_dir} to https://github.com/aertslab/create_cisTarget_databases
create_cistarget_databases_dir='/mnt/beegfs/home/alyne/persistent/create_cisTarget_databases'


#### Score the motifs - This will generate the scores database we will use later on for DEM
${create_cistarget_databases_dir}/create_cistarget_motif_databases.py \
-f ${region_fasta} \
-M ${path_to_motif_collection} \
-m ${motif_list} \
-o ${database_suffix} \
-t ${n_cpu} \
-l \
-s 555


#### Create rankings
${create_cistarget_databases_dir}/convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py \
-i ${database_suffix}.motifs_vs_regions.scores.feather \
-s 555




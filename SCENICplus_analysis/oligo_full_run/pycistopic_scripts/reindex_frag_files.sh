#!/bin/bash
#SBATCH --job-name=ct_reindex
#SBATCH --ntasks=1
#SBATCH --mem=50gb
#SBATCH --array=[1-33]
#SBATCH --time=01:00:00
#SBATCH --partition=dev
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=END
#SBATCH --mail-user=anne-marie.lyne@curie.fr


#initialise and load conda environment
module load utility
source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/pycistopic_polars_3_10_env


#script to do QC for each SAMPLE in parallel
CANCER="oligo"

# cd /mnt/beegfs/home/alyne/stagein/scenicplus_${CANCER}

DATA_DIR=/mnt/beegfs/home/alyne/stagein/scenicplus_${CANCER}

SAMPLE=$(awk -v l="$SLURM_ARRAY_TASK_ID" 'NR==l' /mnt/beegfs/home/alyne/stagein/scenicplus_${CANCER}/sample_list.txt)

cd "$DATA_DIR"/analysis_"$SAMPLE"/CR_output

rm atac_fragments.tsv.gz.tbi

tabix -s 1 -b 2 -e 3 atac_fragments.tsv.gz
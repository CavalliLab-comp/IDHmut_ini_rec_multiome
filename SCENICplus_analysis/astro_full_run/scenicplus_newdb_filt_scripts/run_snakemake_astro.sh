#!/bin/bash
#SBATCH --job-name=scenicplus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=600G
#SBATCH --time=96:00:00
#SBATCH --partition=dev
#SBATCH --output=scplus_%j.out
#SBATCH --error=scplus_%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=anne-marie.lyne@curie.fr


#initialise and load conda environment
module load utility
source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/scenicplus_env

cancer="astro"

cd /mnt/beegfs/home/alyne/stageout/scenicplus_op_newdb_filt_"$cancer"/scplus_pipeline/Snakemake

snakemake --unlock
snakemake --cores 1 --rerun-incomplete --default-resources "tmpdir='/mnt/beegfs/home/alyne/stageout/scenicplus_op_newdb_filt_astro/tmp'"

# snakemake --cores 4 --default-resources "tmpdir='/mnt/beegfs/home/alyne/stageout/scenicplus_op_newdb_filt_astro/tmp'"




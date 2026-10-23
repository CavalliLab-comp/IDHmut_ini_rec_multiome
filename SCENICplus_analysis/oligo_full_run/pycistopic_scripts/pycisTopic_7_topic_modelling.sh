#!/bin/bash
#SBATCH --job-name=ct_top_mod
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=1000gb
#SBATCH --time=72:00:00
#SBATCH --partition=dev
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err
#SBATCH --mail-type=END
#SBATCH --mail-user=anne-marie.lyne@curie.fr


#initialise and load conda environment
module load utility
source /mnt/beegfs/common/apps/miniforge/25.3.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs/home/alyne/persistent/pycistopic_polars_3_10_env

#ensure mallet is downloaded
# cd /mnt/beegfs/home/alyne/persistent
# wget https://github.com/mimno/Mallet/releases/download/v202108/Mallet-202108-bin.tar.gz
# tar -xf Mallet-202108-bin.tar.gz 



CANCER="oligo"
BASE_DIR="/mnt/beegfs/home/alyne"

#create corpus
# mallet_path="Mallet-202108/bin/mallet"
outfile="${BASE_DIR}/stageout/pycisTopic_output_${CANCER}/corpus"
mallet_path="${BASE_DIR}/persistent/Mallet-202108/bin/mallet"
binary_matrix_path="${BASE_DIR}/stageout/pycisTopic_output_${CANCER}/ct_output/binary_accessibility_matrix.mtx"


pycistopic topic_modeling mallet create_corpus -i $binary_matrix_path -o $outfile -m 950 -b $mallet_path



# #run LDA with different numbers of topics

output_prefix="${BASE_DIR}/stageout/pycisTopic_output_${CANCER}/models/tm"
mkdir "${BASE_DIR}/stageout/pycisTopic_output_${CANCER}/models"
mkdir "${BASE_DIR}/stageout/pycisTopic_output_${CANCER}/models/tm"

## run topic modeling
pycistopic topic_modeling mallet run -i $outfile -o $output_prefix -t 30 40 50 60 70 80 \
-p 9 -n 150 -a 50 -A True -e 0.1 -E False -s 555 -m 950 -b $mallet_path -v



regions_path="${BASE_DIR}/stageout/pycisTopic_output_${CANCER}/ct_output/region_names.txt"
cells_path="${BASE_DIR}/stageout/pycisTopic_output_${CANCER}/ct_output/cell_barcodes.txt"
## assemble topic models
pycistopic topic_modeling mallet stats -i $binary_matrix_path -c $cells_path -r $regions_path -o $output_prefix -t 30 40 50 60 70 80 -v





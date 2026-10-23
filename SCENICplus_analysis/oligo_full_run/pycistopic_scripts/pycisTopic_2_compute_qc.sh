#!/bin/bash
#SBATCH --job-name=ct_qc
#SBATCH --ntasks=1
#SBATCH --mem=150gb
#SBATCH --array=[1-33]
#SBATCH --time=02:00:00
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

#make directory for QC output
mkdir /mnt/beegfs/home/alyne/stageout/pycisTopic_output_${CANCER}/qc

#make file with pycistopic qc commands for each SAMPLE
# Nsamp=$(wc -l input_data/SAMPLE_list.txt | awk '{print $1}')

# for ((i=4;i<="$Nsamp";i++)); do

SAMPLE=$(awk -v l="$SLURM_ARRAY_TASK_ID" 'NR==l' /mnt/beegfs/home/alyne/stagein/scenicplus_${CANCER}/sample_list.txt)

#count the number of reads in the fragment file
NREAD=$(zcat "$DATA_DIR"/analysis_${SAMPLE}/CR_output/atac_fragments.tsv.gz | wc -l)

if [ "$NREAD" -ge 120000000 ]; then

#subsample reads to reduce memory
SS_TMP_DIR=/mnt/beegfs/tmp/subsampled_frag_files

#mkdir
mkdir "$SS_TMP_DIR"
mkdir "$SS_TMP_DIR"/analysis_${SAMPLE}
mkdir "$SS_TMP_DIR"/analysis_${SAMPLE}/CR_output

cd "$SS_TMP_DIR"/analysis_${SAMPLE}/CR_output

#subSAMPLE reads to reduce run time
zcat "$DATA_DIR"/analysis_${SAMPLE}/CR_output/atac_fragments.tsv.gz | grep "^#" > headers.txt

N=100000000   # number of lines you want

#sample and sort N random indices
zcat "$DATA_DIR"/analysis_${SAMPLE}/CR_output/atac_fragments.tsv.gz | grep -v "^#" | awk -v k="$N" '
    NR <= k { a[NR] = NR; next }

    NR > k {
        j = int(rand() * NR) + 1
        if (j <= k) a[j] = NR
    }

    END {
        for (i = 1; i <= k; i++) print a[i]
    }
' | sort -n > idx_sorted.txt

#extract these lines
awk '
    NR==FNR { idx[++n] = $1; next }    # read sorted indices (first file)
    BEGIN { i = 1 }                    # ensure i starts at 1 for second file
    {
        if (i > n) exit                # done early when all indices consumed
        if (FNR == idx[i]) {
            print
            i++
        }
    }
' idx_sorted.txt <(zcat "$DATA_DIR"/analysis_${SAMPLE}/CR_output/atac_fragments.tsv.gz | grep -v "^#") > body.txt

#combine files and bgzip
cat headers.txt body.txt | bgzip > atac_fragments.tsv.gz

#remove text files
rm headers.txt
rm body.txt
rm idx_sorted.txt

#index
tabix -s 1 -b 2 -e 3 atac_fragments.tsv.gz

#run pycis topic QC
#NB, requires 150gb per sample
pycistopic qc run -f "$SS_TMP_DIR"/analysis_${SAMPLE}/CR_output/atac_fragments.tsv.gz \
   -r /mnt/beegfs/home/alyne/stageout/pycisTopic_output_${CANCER}/consensus_peak_calling/consensus_regions.bed \
   -t /mnt/beegfs/home/alyne/stagein/required_files/tss.bed \
   -o /mnt/beegfs/home/alyne/stageout/pycisTopic_output_${CANCER}/qc/${SAMPLE} \
   -e polars

echo "finished QC analysis of SAMPLE ${SAMPLE}"

#remove temporary folder
rm -r "$SS_TMP_DIR"/analysis_${SAMPLE}

else

#run on full file
pycistopic qc run -f "$DATA_DIR"/analysis_${SAMPLE}/CR_output/atac_fragments.tsv.gz \
   -r /mnt/beegfs/home/alyne/stageout/pycisTopic_output_${CANCER}/consensus_peak_calling/consensus_regions.bed \
   -t /mnt/beegfs/home/alyne/stagein/required_files/tss.bed \
   -o /mnt/beegfs/home/alyne/stageout/pycisTopic_output_${CANCER}/qc/${SAMPLE} \
   -e polars


fi



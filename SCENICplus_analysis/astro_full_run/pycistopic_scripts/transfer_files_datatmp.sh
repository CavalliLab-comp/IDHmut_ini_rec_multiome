#iterate through each sample and move required files to one folder in /data/tmp

DATA_DIR="/data/users/alyne/gliomult/svn/analyse/results"

DIR_TRANS="/data/tmp/scenicplus_astro"
mkdir $DIR_TRANS
cd $DIR_TRANS

for (( i=1 ; i<=20 ; i++ )); do

SAMP=$(awk -v l=$i 'NR==l' "$DATA_DIR/astro_samples_all.txt")
echo $SAMP

mkdir analysis_"$SAMP"
mkdir analysis_"$SAMP"/CR_output

cp $DATA_DIR/analysis_"$SAMP"/CR_output/atac_fragments.tsv.g* analysis_"$SAMP"/CR_output
cp $DATA_DIR/analysis_"$SAMP"/oligo_filt_final_"$SAMP".rds analysis_"$SAMP"

done




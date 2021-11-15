#!/bin/bash -e
#SBATCH -t 8:00:00 -N 1
#SBATCH --mem=48G
#SBATCH -p blekhman
#SBATCH --mail-type=FAIL

# PROCESSING with DADA2

echo "${PROJECT},dada,start" >> /home/blekhman/shared/compendium/activity.csv

module load singularity

filecount1=$(ls -l /scratch.global/rabdill/bulk/${PROJECT}/fastq/*_1.fastq | wc -l)
filecount2=$(ls -l /scratch.global/rabdill/bulk/${PROJECT}/fastq/*_2.fastq | wc -l)

if [[ $filecount1 -ne $filecount2 ]]
  then
    echo "Running DADA2 in SINGLE-ENDED MODE"
    singularity exec --bind /scratch.global/rabdill/bulk/${PROJECT}:/mnt --bind /home/blekhman/shared/compendium/code:/code docker://blekhmanlab/dada2:1.14.0a Rscript /code/process_forwards.R trim
  else
    echo "Running DADA2 in paired-end mode"
    singularity exec --bind /scratch.global/rabdill/bulk/${PROJECT}:/mnt --bind /home/blekhman/shared/compendium/code:/code docker://blekhmanlab/dada2:1.14.0a Rscript /code/process_paired_end.R trim
fi

echo "Moving results"
cp -R /scratch.global/rabdill/bulk/${PROJECT}/results /home/blekhman/shared/compendium/results/${PROJECT}
echo "DONE"

echo "${PROJECT},dada,end" >> /home/blekhman/shared/compendium/activity.csv

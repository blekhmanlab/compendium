#!/bin/bash -e
#SBATCH -t 16:00:00 -N 1
#SBATCH --mem=32G
#SBATCH --ntasks=8
#SBATCH -p blekhman
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rabdill@umn.edu

# PROCESSING with DADA2

echo "${PROJECT},matchdada,start" >> /home/blekhman/shared/compendium/activity.csv

module load singularity

echo "Running DADA2 in paired-end mode"
singularity exec --bind /scratch.global/rabdill/bulk/${PROJECT}:/mnt --bind /home/blekhman/shared/compendium/code:/code /home/blekhman/shared/compendium/code/dada2_1.14.0a.sif Rscript /code/process_paired_end.R matchids

echo "Moving results"
cp -R /scratch.global/rabdill/bulk/${PROJECT}/results /home/blekhman/shared/compendium/results/${PROJECT}
echo "DONE"

echo "${PROJECT},matchdada,end" >> /home/blekhman/shared/compendium/activity.csv

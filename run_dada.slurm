#!/bin/bash -e
#SBATCH -t 3000 -N 1
#SBATCH --mem=48G
#SBATCH -p blekhman
#SBATCH --mail-type=FAIL

# PROCESSING with DADA2

echo "${PROJECT},dada,start" >> /home/blekhman/shared/compendium/activity.csv

module load singularity

cmd=$(cat /scratch.global/rabdill/bulk/${PROJECT}/learnt/shi7_cmd.sh)

if [[ "$cmd" == *"-SE"* ]]
  then
    echo "Running DADA2 in SINGLE-ENDED MODE"
    singularity exec --bind /scratch.global/rabdill/bulk/${PROJECT}:/mnt docker://blekhmanlab/dada2:1.14.0a Rscript /home/blekhman/shared/compendium/code/process_forwards.R
  else
    echo "Running DADA2 in paired-end mode"
    singularity exec --bind /scratch.global/rabdill/bulk/${PROJECT}:/mnt docker://blekhmanlab/dada2:1.14.0a Rscript /home/blekhman/shared/compendium/code/process_paired_end.R
fi

echo "Moving results"
cp -R /scratch.global/rabdill/bulk/${PROJECT}/results /home/blekhman/shared/compendium/results/${PROJECT}
echo "DONE"

echo "${PROJECT},dada,end" >> /home/blekhman/shared/compendium/activity.csv

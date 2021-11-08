#!/bin/bash -e
#SBATCH -t 24:00:00 -N 1
#SBATCH --mem=2G
#SBATCH -p blekhman
#SBATCH --mail-type=FAIL

# This script expects a directory matching the project name
# to be present at /home/blekhman/shared/compendium/accession_lists/
# with a file in it called "SraAccList.txt" listing all the samples
# to be downloaded from SRA.

if [ -z "${PROJECT}" ]
then
      echo "PROJECT variable is EMPTY. Bailing"
      exit 1
fi

cd /scratch.global/rabdill/bulk/${PROJECT}

module load sratoolkit # currently 2.8.2

echo "Converting files to FASTQ"
for sample in $(cat SraAccList.txt)
do
    echo "On sample: $sample"
    fastq-dump -v -I --split-files $sample -O /scratch.global/rabdill/bulk/${PROJECT}/fastq
done

# create a "samples" file listing ONLY the samples that were successfully processed
cd fastq
ls *_1.fastq | cut -f1 -d "_" > ../samples.txt
echo "DONE. Submitting DADA2 job."
echo "${PROJECT},download,end" >> /home/blekhman/shared/compendium/activity.csv

cd /home/blekhman/shared/compendium/logs
sbatch --export=PROJECT=${PROJECT} --job-name=${PROJECT}dada -o ${PROJECT}dada.log ../code/run_dada.slurm

#!/bin/bash -e
#SBATCH -t 32:00:00 -N 1
#SBATCH --ntasks=4
#SBATCH --mem=1G
#SBATCH -p blekhman
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rabdill@umn.edu

# This script expects a directory matching the project name
# to be present at /home/blekhman/shared/compendium/accession_lists/
# with a file in it called "SraAccList.txt" listing all the samples
# to be downloaded from SRA.

if [ -z "${PROJECT}" ]
then
      echo "PROJECT variable is EMPTY. Bailing"
      exit 1
fi

echo "${PROJECT},download,start" >> /home/blekhman/shared/compendium/activity.csv

cp -R /home/blekhman/shared/compendium/accession_lists/${PROJECT} /scratch.global/rabdill/bulk/${PROJECT}
cd /scratch.global/rabdill/bulk/${PROJECT}
mkdir -p temp
mkdir -p results
mkdir -p intermediate

echo "Fetching SRA files"
num=1
for sample in $(cat SraAccList.txt)
do
    echo "On sample ${num}: $sample"
    if [[ -f "/scratch.global/rabdill/bulk/${PROJECT}/fastq/${sample}_1.fastq" ]]; then
        echo "Already found! Continuing."
    else
        /home/blekhman/shared/compendium/code/sratoolkit.2.11.3-centos_linux64/bin/fasterq-dump $sample -e 4 -S -O /scratch.global/rabdill/bulk/${PROJECT}/fastq || echo "SOMETHING WENT WRONG HERE. Continuing..."
    fi
    num=$(($num + 1))
done

echo "DONE. Submitting dada job."

echo "${PROJECT},download,end" >> /home/blekhman/shared/compendium/activity.csv

cd fastq
ls *_1.fastq | cut -f1 -d "_" > ../samples.txt

echo "${PROJECT},unpack,start" >> /home/blekhman/shared/compendium/activity.csv
echo "${PROJECT},unpack,end" >> /home/blekhman/shared/compendium/activity.csv


cd /home/blekhman/shared/compendium/logs
sbatch --export=PROJECT=${PROJECT} --job-name=${PROJECT}dada -o ${PROJECT}dada.log ../code/run_dada.slurm

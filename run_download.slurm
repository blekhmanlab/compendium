#!/bin/bash -e
#SBATCH -t 12:00:00 -N 1
#SBATCH --mem=2G
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
mkdir fastq
mkdir temp
mkdir results
mkdir intermediate

#module load aspera
module load sratoolkit # currently 2.8.2

# configuration: vdb-config -i

# find a list of samples, pick one, then go to
# the PROJECT page for it. Click on "SRA Experiments"
# and then export all the accession numbers to a file
# called SraAccList.txt

echo "Fetching SRA files"
prefetch --option-file SraAccList.txt

echo "DONE. Submitting unpacking job."
echo "${PROJECT},download,end" >> /home/blekhman/shared/compendium/activity.csv

cd /home/blekhman/shared/compendium/logs
sbatch --export=PROJECT=$PROJECT --job-name=${PROJECT}unpack -o ${PROJECT}unpack.log ../code/download_unpack.slurm

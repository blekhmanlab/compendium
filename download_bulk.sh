#!/bin/sh -e

# This script expects a directory matching the project name
# to be present at /scratch.global/rabdill/bulk/
# with a file in it called SraAccList.txt listing all the samples
# to be downloaded from SRA.

mv ~/shithouse/accession_lists/${1} /scratch.global/rabdill/bulk/
cd /scratch.global/rabdill/bulk/${1}
mkdir fastq
mkdir intermediate

#module load aspera
module load sratoolkit # currently 2.8.2

# configuration: vdb-config -i

# find a list of samples, pick one, then go to
# the PROJECT page for it. Click on "SRA Experiments"
# and then export all the accession numbers to a file
# called SraAccList.txt

prefetch --option-file SraAccList.txt

for sample in $(cat SraAccList.txt)
do
    echo "On sample: $sample"
    fastq-dump -I --split-files $sample -O /scratch.global/rabdill/bulk/${1}/fastq
done

# create a "samples" file listing ONLY the samples that were successfully processed
cd fastq
ls *_1.fastq | cut -f1 -d "_" > ../samples.txt
echo "DONE"

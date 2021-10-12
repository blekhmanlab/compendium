#!/bin/bash -e
#SBATCH -t 960 -N 1
#SBATCH --mem=2G
#SBATCH -p blekhman
#SBATCH --mail-type=FAIL

echo "${PROJECT},trim,start" >> ~/shithouse/activity.csv

mkdir /scratch.global/rabdill/bulk/${PROJECT}/results
mkdir /scratch.global/rabdill/bulk/${PROJECT}/trimmed # for shi7
mkdir /scratch.global/rabdill/bulk/${PROJECT}/intermediate # for dada2
mkdir /scratch.global/rabdill/bulk/${PROJECT}/temp

# make sure there isn't a readlength record from some previous failed run
rm -f /scratch.global/rabdill/bulk/${PROJECT}/temp/readlength.raw.txt

# record read lengths
echo "Assessing read length..."
for filename in /scratch.global/rabdill/bulk/${PROJECT}/fastq/*.fastq; do
    [ -e "$filename" ] || continue # make sure something matches
    cat $filename | perl -ne '$s=<>;<>;<>;chomp($s);print length($s)."\n";' >> /scratch.global/rabdill/bulk/${PROJECT}/temp/readlength.raw.txt
done
cat /scratch.global/rabdill/bulk/${PROJECT}/temp/readlength.raw.txt | uniq -c >> /scratch.global/rabdill/bulk/${PROJECT}/results/readlength.txt

# TRIMMING with shi7
# (all this weird formatting is because we have to run shi7 SPECIFICALLY
#  with python 2, rather than the python 3 that the system defaults to,
#  because of an outstanding bug)
cd ~/shithouse/shi7_0.9.9/
echo "Learning shi7 parameters..."
python2 shi7_learning.py -i /scratch.global/rabdill/bulk/${PROJECT}/fastq -o /scratch.global/rabdill/bulk/${PROJECT}/learnt
echo "Running shi7..."
a=$(cat /scratch.global/rabdill/bulk/${PROJECT}/learnt/shi7_cmd.sh)

$(echo python2 $a --flash False --convert_fasta False --combine_fasta False --output /scratch.global/rabdill/bulk/${PROJECT}/trimmed)

echo "DONE."
echo "${PROJECT},trim,end" >> ~/shithouse/activity.csv

files=$(ls -l /scratch.global/rabdill/bulk/${PROJECT}/trimmed | wc -l)
mem=16
case $files in
[1-9][0-9]) # 99
    mem=16
    ;;
[1-2][0-9][0-9]) # 299
    mem=48
    ;;
[3-4][0-9][0-9]) # 399
    mem=64
    ;;
[4-6][0-9][0-9]) # 699
    mem=80
    ;;
[7-9][0-9][0-9]) # 999
    mem=96
    ;;
1[0-9][0-9][0-9]) # 1999
    mem=112
    ;;
2[0-4][0-9][0-9]) # 2499
    mem=128
    ;;
2[5-9][0-9][0-9]) # 2999
    mem=144
    ;;
3[0-4][0-9][0-9]) # 3499
    mem=160
    ;;
3[5-9][0-9][0-9])  # 3999
    mem=176
    ;;
4[0-9][0-9][0-9]) # 4999
    mem=192
    ;;
[1-9]*)
    mem=256
    ;;
*)
    echo "Unable to determine memory requirements for count '$files'; something's wrong."
    exit
    ;;
esac

echo "Starting DADA2 run with ${mem} GB of memory."
cd ~/shithouse/logs
sbatch --mem=${mem}G --export=PROJECT=${PROJECT} --job-name=${PROJECT}dada -o ${PROJECT}dada.log ../code/run_dada.slurm

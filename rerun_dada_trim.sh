#!/bin/sh

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

# Removing reverse reads
bash remove_reverse_reads.sh $1

rm -rf /home/blekhman/shared/compendium/results/${1}
rm /scratch.global/rabdill/bulk/${1}/results/*
rm /scratch.global/rabdill/bulk/${1}/temp/*
rm /scratch.global/rabdill/bulk/${1}/intermediate/*

cd /home/blekhman/shared/compendium/logs
sbatch --export=PROJECT=$1 --job-name=${1}trimdada -o ${1}trimdada.log ../code/run_trim_dada.slurm

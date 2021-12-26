#!/bin/sh

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

rm -rf /home/blekhman/shared/compendium/results/${1}
rm /scratch.global/rabdill/bulk/${1}/results/*
rm /scratch.global/rabdill/bulk/${1}/temp/*
rm /scratch.global/rabdill/bulk/${1}/intermediate/*

cd /home/blekhman/shared/compendium/logs
sbatch --export=PROJECT=$1 --job-name=${1}dadamatch -o ${1}matchdada.log ../code/run_match_dada.slurm

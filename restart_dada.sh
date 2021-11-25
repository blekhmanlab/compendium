#!/bin/sh

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

rm /scratch.global/rabdill/bulk/${1}/results/*
rm /scratch.global/rabdill/bulk/${1}/temp/*
rm /scratch.global/rabdill/bulk/${1}/intermediate/*

cd /home/blekhman/shared/compendium/logs
sbatch --export=PROJECT=$1 --job-name=${1}dada -o ${1}dada.log ../code/run_dada.slurm

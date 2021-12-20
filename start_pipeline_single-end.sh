#!/bin/sh

cd /home/blekhman/shared/compendium/logs
sbatch --export=PROJECT=$1 --job-name=${1}download -o ${1}download.log ../code/run_fasterq_single-end.slurm

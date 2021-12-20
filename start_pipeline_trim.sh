#!/bin/sh

# This one starts the pipeline and jumps straight to
# running the "trim" version of DADA2
cd /home/blekhman/shared/compendium/logs
sbatch --export=PROJECT=$1 --job-name=${1}download -o ${1}download.log ../code/run_fasterq_trim.slurm

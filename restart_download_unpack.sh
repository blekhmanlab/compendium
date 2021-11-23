#!/bin/bash

cd /home/blekhman/shared/compendium/logs
sbatch --export=PROJECT=$1 --job-name=${1}download -o ${1}download.log ../code/download_unpack.slurm

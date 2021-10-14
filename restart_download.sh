#!/bin/sh
# NOTE: This should be run from within the logs directory.
# It's probably best to call this from inside the
# "run_restart_download.slurm" file.
sbatch --export=PROJECT=$1 --job-name=${1}redownload -o ${1}redownload.log ../code/run_restart_download.slurm

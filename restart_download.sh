#!/bin/sh
# NOTE: This should be run from within the logs directory.
cd /home/blekhman/shared/compendium/logs

sbatch --export=PROJECT=$1 --job-name=${1}redownload -o ${1}redownload.log ../code/run_restart_download.slurm

#!/bin/sh

sbatch --export=PROJECT=$1 --job-name=${1}redownload -o ${1}redownload.log ../code/run_restart_download.slurm

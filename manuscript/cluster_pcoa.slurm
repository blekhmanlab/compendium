#!/bin/bash -e
#SBATCH -t 48:00:00 -N 1
#SBATCH --ntasks=4
#SBATCH --mem=240G
#SBATCH -p blekhman
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rabdill@umn.edu

module load R/4.0.4
Rscript cluster_pcoa.R


# sbatch --job-name=PCOA -o pcoa.log cluster_pcoa.slurm
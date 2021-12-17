#!/bin/bash -e
#SBATCH -t 12:00:00 -N 1
#SBATCH --mem=16G
#SBATCH -p blekhman
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rabdill@umn.edu

# PROCESSING with DADA2

cd /home/blekhman/shared/compendium/code
python3 consolidate_taxa.py
# python3 combine_studies.py
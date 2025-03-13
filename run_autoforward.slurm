#!/bin/bash -e
#SBATCH -t 00:10:00 -N 1
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH -A pi-blekhman
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rabdill@uchicago.edu

module load python
source venv/bin/activate
python main.py autoforward

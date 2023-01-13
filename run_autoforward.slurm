#!/bin/bash -e
#SBATCH -t 00:10:00 -N 1
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH -p blekhman
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rabdill@uchicago.edu

module load python

if [ ! -d "./venv" ]; then
    python -m venv venv
fi

source venv/bin/activate

if [ ! -f "./venv/bin/pytest" ]; then
    pip install --upgrade pip
    pip install -r requirements.txt
fi

python main.py autoforward

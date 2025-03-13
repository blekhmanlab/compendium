#!/bin/bash -e

# project ID
gotID=$(echo -n "$1" | wc -c) # count how long this string is
if [ $gotID -lt 5 ]; then
    echo "ProjectID parameter not long enough"
    exit 1
fi

cd "$1"

if [ ! -d "./venv" ]; then
    echo "Setting up virtalenv"
    python -m venv venv
fi

source venv/bin/activate

if [ ! -f "./venv/bin/snakemake" ]; then
    echo "Running pip"
    pip install --upgrade pip
    pip install -r requirements.txt
fi

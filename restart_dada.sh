#!/bin/sh

if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters: need project, hours and memory"
    exit 1
fi

proj=$1
hours=$2
gb=$3

# make sure we got a valid NEW value
re='^[0-9]+$'
if ! [[ $hours =~ $re ]] ; then
    echo "error: New duration not a number"
    exit 1
fi

if ! [[ $gb =~ $re ]] ; then
    echo "error: New memory not a number"
    exit 1
fi

echo "Re-running DADA with ${hours} hours and ${gb} GB of memory"

filename=$(date +%s).slurm # timestamp is filename

# check valid filename
if [[ ${#filename} -lt 8 ]]; then
    echo "Automatically generated filename not valid. Bailing."
    exit 1
fi

echo "#!/bin/bash -e" > $filename
echo "#SBATCH -t ${hours}:00:00 -N 1" >> $filename
echo "#SBATCH --mem=${gb}G" >> $filename
cat dada.template >> $filename

rm /scratch.global/rabdill/bulk/${1}/results/*
rm /scratch.global/rabdill/bulk/${1}/temp/*
rm /scratch.global/rabdill/bulk/${1}/intermediate/*

cd /home/blekhman/shared/compendium/logs
sbatch --export=PROJECT=$1 --job-name=${1}dada -o ${1}dada.log ../code/${filename}

rm ../code/$filename

#!/bin/sh

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi
if [ ${#1} -lt 4 ]; then
    echo "Project name isn't long enough. Bailing"
    exit 1
fi

rm -rf /scratch.global/rabdill/bulk/${1}
rm -rf /home/blekhman/shared/compendium/results/${1}
rm /home/blekhman/shared/compendium/results/taxa_files/${1}_consolidated.tsv

cat ../activity.csv | grep -v ${1}, >> ../activity.csv
cat to_ignore.csv | grep -v -w ${1} >> to_ignore.csv

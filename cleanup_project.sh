#!/bin/bash -e

head /home/blekhman/shared/compendium/results/${1}/ASVs_taxonomy.tsv
head /home/blekhman/shared/compendium/results/${1}/summary.tsv
echo "\\nClean up project ${1}? Type one of the following options."
echo "\\n\\ny = 'Yes, we want to include this, clean up the extra files.'"
echo "n = 'Nevermind, leave everything where it is.'"
echo "DISCARD = 'Something is wrong with this, GET RID OF IT PERMANENTLY.'\\n"
read doit

if [[ "$doit" == "DISCARD" ]]; then
    echo 'Are you SURE you want to delete this project? (y/n)'
    read doit
    if [[ "$doit" != "y" ]]; then
        echo "Ok, not doing anything."
        exit 0
    fi
    if [[ ${#1} -lt 4 ]]; then
        echo "Length of project ID is less than 4. Skipping to avoid doing something bad."
        exit 1
    fi

    echo 'DISCARDING.'
    echo "/scratch.global/rabdill/bulk/${1}"
    #rm -rf /scratch.global/rabdill/bulk/${1}
    echo ${1} >> /home/blekhman/shared/compendium/code/to_ignore.csv
    exit 0
fi

if [[ "$doit" != "y" ]]; then
    echo 'Skipping.'
    exit 0
fi

echo "${1},archive,start" >> /home/blekhman/shared/compendium/activity.csv

cd /scratch.global/rabdill/bulk/${1}
rm -rf temp
rm -rf intermediate
rm -rf fastq
rm -rf results
cd ..
tar -zcvf ${1}.tar.gz ${1}
mv ${1}.tar.gz /home/blekhman/shared/compendium/archived/
rm -rf ${1}
echo "Done"
echo "${1},archive,end" >> /home/blekhman/shared/compendium/activity.csv

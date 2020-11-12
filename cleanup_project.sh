#!/bin/bash -e

head ~/shithouse/results/${1}/ASVs_taxonomy.tsv
head ~/shithouse/results/${1}/summary.tsv
echo "Clean up project ${1}?"
read doit
if [[ $doit!='y' ]]; then
    echo 'Skipping.'
    exit 0
fi

echo "${1},archive,start" >> ~/shithouse/activity.csv

cd /scratch.global/rabdill/bulk/${1}
rm -rf temp
mv trimmed/shi7.log .
mv learnt/shi7_learning.log .
mv learnt/shi7_cmd.sh .
rm -rf trimmed
rm -rf intermediate
rm -rf learnt
rm -rf fastq
rm -rf results
cd ..
tar -zcvf ${1}.tar.gz ${1}
mv ${1}.tar.gz ~/shithouse/archived/
rm -rf ${1}
echo "Done"
echo "${1},archive,end" >> ~/shithouse/activity.csv

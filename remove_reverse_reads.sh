#!/bin/sh

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi
if [ ${#1} -lt 4 ]; then
    echo "Project name isn't long enough. Bailing"
    exit 1
fi

rm -f /scratch.global/rabdill/bulk/${1}/fastq/*_2.fastq

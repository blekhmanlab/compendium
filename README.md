# shithouse
**Large-scale characterization of the human microbiome using machine learning**

```sh
cd ~/shithouse/logs
# downloads fastq files:
qsub ../code/download_bulk.pbs -v PROJECT=PRJNA493625

# extract ASV table for study
qsub ../code/process_project.pbs -v PROJECT=PRJNA493625

# clean up the extra files
bash ../code/cleanup_project.sh PRJNA493625
```

This should (ideally) result in a file at `~/shithouse/results/summaries/PRJNA493625_summary.tsv` and an ASV table at `~/shithouse/results/PRJNA493625.tsv`

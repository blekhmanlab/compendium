# shithouse
**Large-scale characterization of the human microbiome using machine learning**

```sh
cd ~/shithouse
# downloads fastq files:
qsub download_bulk.pbs -v PROJECT=PRJNA493625

# extract ASV table for study
qsub process_project.pbs -v PROJECT=PRJNA493625
```

This should (ideally) result in a file at `~/shithouse/results/summaries/PRJNA493625_summary.tsv` and an ASV table at `~/shithouse/results/PRJNA493625.tsv`

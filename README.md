# shithouse
**Large-scale characterization of the human microbiome using machine learning**

```sh
cd ~/shithouse/logs
# downloads fastq files:
qsub ../code/download_bulk.pbs -v PROJECT=PRJNA513489

# trim reads (the download script should submit this automatically now)
#qsub ../code/run_shi7.pbs -v PROJECT=PRJNA493625

# extract ASV table for study
qsub ../code/run_dada.pbs -v PROJECT=PRJNA666641

# clean up the extra files
bash ../code/cleanup_project.sh PRJNA493625
```

This should (ideally) result in a file at `~/shithouse/results/summaries/PRJNA493625_summary.tsv` and an ASV table at `~/shithouse/results/PRJNA493625.tsv`

## Troubleshooting
If DADA2 fails with
```
Error in filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads,  : 
  Some input files do not exist.
```
Check on the file at `$PROJECT/learnt/shi7_cmd.sh`. If it was run with the `-SE` flag, the data probably isn't paired-end but should be looked at more closely.
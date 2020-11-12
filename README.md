# shithouse
**Large-scale characterization of the human microbiome using machine learning**

Kicking off the download script should start a chain of jobs that process the whole study:
```sh
cd ~/shithouse/logs
qsub ../code/download_bulk.pbs -v PROJECT=PRJNA547558
```

Then eyeball the results to make sure everything finished properly:

```sh
tail PRJNA513489dada.*
head ../results/PRJNA513489/summary.tsv
head ../results/PRJNA513489/ASVs_taxonomy.tsv
```

Then clean up the extra files:

```sh
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
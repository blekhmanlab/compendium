# shithouse
**Large-scale characterization of the human microbiome using machine learning**

## Downloading metadata from NCBI

The `main.py` file has multiple functions that accommodate different steps in the process of pulling data out of the Sequence Read Archive and putting it into the shithouse database.

* **`load_xml()`:** This loads an XML file **exported from a BioSample search** and puts the sample metadata into the database.
* **`find_runs()`:** This step uses the sample metadata (*already in the database*) to find *runs* ("SRR" codes), which are the entities actually associated with sequencing data.
* **`write_lists()`:** This function fetches the SRA projects and generates text files containing runs for each project. **These accession lists are the input for the processing pipeline.**

## Processing SRA files
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

You can check on projects that got stalled in the pipeline by running:
```sh
python3 ../code/review_activites.py
```

## Troubleshooting
If DADA2 fails with
```
Error in filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads,  :
  Some input files do not exist.
```
Check on the file at `$PROJECT/learnt/shi7_cmd.sh`. If it was run with the `-SE` flag, the data probably isn't paired-end but should be looked at more closely.

# Human Microbiome Compendium
**Large-scale characterization of the human microbiome using machine learning**

## Downloading metadata from NCBI

The `main.py` file has multiple functions that accommodate different steps in the process of pulling data out of the Sequence Read Archive and putting it into the samples database.

* **`load_xml()`:** This loads an XML file **exported from a BioSample search** and puts the sample metadata into the database.
* **`find_runs()`:** This step uses the sample metadata (*already in the database*) to find *runs* ("SRR" codes), which are the entities actually associated with sequencing data.
* **`write_lists()`:** This function fetches the SRA projects and generates text files containing runs for each project. **These accession lists are the input for the processing pipeline.**

## Processing SRA files
Kicking off the download script should start a chain of jobs that process the whole study:
```sh
cd /home/blekhman/shared/compendium/code
bash start_pipeline.sh PRJNA547558
```
When it's done, clean up the extra files:

```sh
bash cleanup_project.sh PRJNA493625
```

You can check on projects that got stalled in the pipeline by running:
```sh
python3 review_activites.py
```

## Troubleshooting

1. If DADA2 fails with
```
Error in filterAndTrim(forward_reads, filtered_forward_reads, reverse_reads,  :
  Some input files do not exist.
```
Check on the file at `$PROJECT/learnt/shi7_cmd.sh`. If it was run with the `-SE` flag, the data probably isn't paired-end but should be looked at more closely.

2. If a high percentage of reads were removed by the "remove chimeras" step, it probably means the primers were still attached to the reads and the trimming step didn't catch them. Flag this and we can return to it later.
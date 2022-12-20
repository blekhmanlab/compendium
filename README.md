# HMC Metadata Ingest Service

## Downloading metadata from NCBI
We currently extract relevant samples from search results on [the BioSample website](https://www.ncbi.nlm.nih.gov/biosample) using this query:
```
txid408170[Organism:noexp] AND ("public"[filter] AND "biosample sra"[filter])
```
Once the results are displayed, select "Send to" > "File" with the "Full XML (text)" format. This will probably take a while.

## Ingesting
The `main.py` file has multiple functions that accommodate different steps in the process of pulling data out of the Sequence Read Archive and putting it into the samples database.

* **`load_xml()`:** This loads an XML file **exported from a BioSample search** and puts the sample metadata into the database.
* **`find_runs()`:** This step uses the sample metadata (*already in the database*) to find *runs* ("SRR" codes), which are the entities actually associated with sequencing data.
* **`write_lists()`:** This function fetches the SRA projects and generates text files containing runs for each project. **These accession lists are the input for the processing pipeline.**

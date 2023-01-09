# HMC Metadata Ingest Service

## Configuration
Create a file called `config.py` with the following values, all of which are passed in API requests to NCBI:

* **Tool**: The name of the application that's collecting data.
* **Email**: A contact address NCBI administrators can reach.
* **Key**: An NCBI API key.

## Project states
The `status` table tracks the state of all projects that have been referred to by the CLI in some way. (Projects that had their metadata loaded, but were never processed, are not there.) The possible conditions are:

* `accession_list_created`: A directory has been created for the project and a list of its samples has been added.
* `running`: A job has been submitted to SLURM to be processed
* `to_re_run`: A paired-end project has been flagged for getting re-run as single-end. (This changes to "running" once it's actually been submitted, which happens almost immediately.)
* `failed`: *Terminal status.* Indicates the project was discarded.
* `complete`: Project is done and had acceptable results that were successfully loaded into the database. **Not** a terminal status.
* `archived`: The project's results were stored in a tarball.
* `done`: *Terminal status.* Indicates the project's results were loaded into the database, its files archived, and its other files cleaned up and deleted.

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

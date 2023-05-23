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

## Commands

The manager application is invoked by running `python main.py` with command-line parameters indicating which operation to perform (e.g. `python main.py runs 2000`). The following options are currently available:

* **xml**: Loads data exported from the BioProject portal (see "Downloading metadata from NCBI" above for details). Data is pulled from this file and  Two required parameters:
  * **taxon** – the NCBI taxon ID used in the search (e.g. txid408170)
  * **filename** - the path to the XML file to be parsed.
* **runs**: Queries the compendium database for samples that have an SRS (sample) number, but not an SRR (run) number. This list is then sent to the NCBI eUtils API to retrieve the runs. The only parameter is a limit on how many samples to evaluate (default 2000).
* **lists**: Generates project-level text files listing SRR numbers for all samples within the project. Each is placed in the `accession_lists` directory and named for the BioProject ID, e.g. `PRJNA12345.txt`. Two optional (positional) parameters:
  * **min_samples**: What is the minimum number of samples a project must have before it's recorded? (default 50)
  * **max_samples**: What is the MAXIMUM number of samples a project must have before it's recorded? (default 100000)
* **runit**: Initializes the processing pipeline for a single project and starts the pipeline. This is generally used to start the pipeline *for the first time*, because it creates a new directory for the project and pulls in all the necessary pipeline code. It will throw a warning if the project has been initialized before, but you can safely proceed if the previous run has been removed and you really do want to start over completely. One parameter:
  * **project**: The BioProject ID of a project for which the metadata is already in our database. (This will retrieve and process the actual files.)
* **discard**: Throws out any computational results from a single project and records in the database that the project should not be re-attempted. Will also prompt you for a brief explanation of why it should be skipped.
  * **project**: The BioProject ID of the project to be thrown out. Required.
* **again**: Submits a new slurm job to restart the snakemake pipeline for a single project. *This command assumes the pipeline has already been configured.* Used mostly for situations in which a project stalled for reasons that have been remediated. One parameter:
  * **project**: The BioProject ID of the project to be restarted. Required.
* **status**: Retrieves the pipeline progress of a single project and prints a report for the user. One parameter:
  * **project**: The BioProject ID of the project to check on. Required.
* **eval**: Checks the progress of a single study. If it's completed the pipeline, it will evaluate the results and prompt the user to confirm that the project should either be saved and finalized, OR should be re-run with different parameters. One parameter:
  * **project**: The BioProject ID of the project to check on. Required.
* **compendium**: Retrieves a summary of progress in processing projects for the compendium and prints a report.
* **summary**: Retrieves a list of all projects currently in progress and prints a report about each.
* **FORWARD**: Iterates through projects with things that need to be addressed and prompts the user to approve the actions. Projects that are still running, or that failed for unknown reasons, simply have their status printed.
* **autoforward**: Similar to FORWARD, but automatically approves actions that need to be taken. If projects are completed, the application will then search for new projects to start.

## Ingesting
The `main.py` file has multiple functions that accommodate different steps in the process of pulling data out of the Sequence Read Archive and putting it into the samples database.

* **`load_xml()`:** This loads an XML file **exported from a BioSample search** and puts the sample metadata into the database.
* **`find_runs()`:** This step uses the sample metadata (*already in the database*) to find *runs* ("SRR" codes), which are the entities actually associated with sequencing data.
* **`write_lists()`:** This function fetches the SRA projects and generates text files containing runs for each project. **These accession lists are the input for the processing pipeline.**

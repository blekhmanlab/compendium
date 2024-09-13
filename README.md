![Human Microbiome Compendium logo](http://blekhmanlab.org/images/compendium.png "Human Microbiome Compendium")

Welcome to the central repository for the Human Microbiome Compendium, an ongoing project to process and integrate hundreds of thousands of human microbiome samples using publicly available data from the members of the [International Nucleotide Sequence Database Collaboration](https://www.insdc.org/). In short, we retrieve sequencing data from the **BioProject and Sequence Read Archive** databases, process it with a uniform pipeline, and combine it into **one large dataset** that can be used as a training set for machine-learning models, a conventional dataset for microbiome analysis, or as additional context for samples of your own.

**If you have feedback, questions or issues with the HMC dataset or code, you're in the right place:**

* For bug fixes, feature requests and suggestions, please [submit a new Issue](https://github.com/blekhmanlab/compendium/issues).
* For troubleshooting help or questions about a particular use case, please [submit a new Discussion topic](https://github.com/blekhmanlab/compendium/discussions).
* For privacy or security concerns, please see our [security policy](/SECURITY.md).

Information about the project is spread across several locations:

* Our website at **[microbiomap.org](https://microbiomap.org)** provides information about the project and the most up-to-date links for announcements, releases and publications.
    * The [`compendium_website`](https://github.com/blekhmanlab/compendium_website) repository contains the code for microbiomap.org.
* This repository contains the code for the compendium management software developed to automate steps for processing and quality control.
* The [`snakemake-compendium`](https://github.com/blekhmanlab/snakemake-compendium) repository contains the pipeline code used to process samples. Compendium Manager launches individual instances of this pipeline for each BioProject.
* The [`MicroBioMap`](https://github.com/seandavi/MicroBioMap) repository contains the code for the R package developed to streamline the retrieval and loading of compendium data.

Please note, the application and pipeline code in these repositories represents the most up-to-date processes in use by the compendium for what will be the 2.0 release. Code used for processing all 1.x releases of the compendium is [archived on Zenodo](https://doi.org/10.5281/zenodo.13733482).

---

# HMC Compendium Manager

This is the command-line utility being developed for use with the Human Microbiome Compendium, but it may be useful to others looking to deal with bulk genomic data from NCBI. It ingests metadata about BioSamples, pulls together enough information to download the deposited FASTQ files from the Sequence Read Archive, and deploys individual Snakemake pipelines to process each project separately.

**This software has not yet had a public 1.0 release,**  nor has it been tested on anything other than our very specific use case—please keep this in mind when using this tool or submitting feedback.

## Structure

This "ingest service" is designed primarily to launch Snakemake pipelines, monitor their progress, and collect their results upon completion. However, there are multiple steps that happen before and after the pipelines are deployed. Broadly, these are the steps that make up the complete workflow:

1. Search results are exported from the BioSample web portal (see below) into an XML file. This file forms the core of the data used by the ingest service, which parses the XML, extracts the relevant data and saves it to a local SQLite database.
1. The application uses the NCBI "E-utilities" to retrieve enough information about each BioSample to associate each one with entries in the Sequence Read Archive. (This can take a long time for large collections of samples.)
1. The application reviews the list of projects and locates one that has not yet been processed. A copy of the Snakemake pipeline is created for this project.
1. The application submits a batch job to the HPC scheduler that launches Snakemake. This job is Snakemake's "supervisor" process—it runs for the duration of the pipeline, and Snakemake uses it to monitor progress through the pipeline. Snakemake submits its own batch jobs for each sample in each step of the pipeline. **The application itself now exits.**
1. The final step of the Snakemake pipeline, if it completes successfully, is to **submit a new batch job** that runs the application's "forward" command. This command prompts the application to check the status of all currently running projects and update its records.
1. If the application observes any projects that have completed since it last checked, the application loads summary information about the results and validates that the project's results are acceptable—whether the proportion of chimeric reads is too high, for example, or whether a suspiciously low proportion of forward reads were matched to reverse reads. If paired reads cannot be reliably merged, reverse reads are discarded and the project is reprocessed as a single-end dataset.
1. If the project passes all quality control checks, the results are parsed out of the pipeline's text file output and loaded into the local database.

## Installation

This application requires Python 3 and has been tested with Python 3.9.6. The CLI should work as expected on MacOS and Linux, but the Snakemake pipelines it invokes for processing the raw data will likely only work on Linux. Please see [the associated pipeline repository](https://github.com/blekhmanlab/snakemake-compendium) for details.

From the command-line of your choice, run these commands:

```sh
git clone git@github.com:blekhmanlab/compendium-manager.git
cd compendium-manager
###### (optional: establish a virtual environment)
python -m venv venv
source venv/bin/activate
######
pip install -r requirements.txt
```

## Configuration
Copy the `config_template.py` file and name it `config.py`. All of the options are editable, but the following variables *must* be set so they can be appended to API requests sent to NCBI:

* **Tool**: The name of the application that's collecting data.
* **Email**: A contact address NCBI administrators can reach.
* **Key**: An [NCBI API key](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/), available for free.

There are many other options that can be tweaked in the config file; comments around each value explain their use.

## Commands

For now, the manager application is invoked by running `python main.py` with command-line parameters indicating which operation to perform (e.g. `python main.py runs 2000`).

### Compendium-level commands
* **`compendium`**: Retrieves a summary of progress in processing projects for the compendium and prints a report.
* **`summary`**: Retrieves a list of all projects currently in progress and prints a report about each.
* **`FORWARD`**: Iterates through projects with things that need to be addressed and prompts the user to approve the actions. Projects that are still running, or that failed for unknown reasons, simply have their status printed.
* **`autoforward`**: Similar to FORWARD, but automatically approves actions that need to be taken. If projects are completed, the application will then search for new projects to start.

### Project-level commands

* **`status`**: Retrieves the pipeline progress of a single project and prints a report for the user. One parameter, the BioProject ID of the project to check on.
* **`eval`**: Checks the progress of a single study. If it's completed the pipeline, it will evaluate the results and prompt the user to confirm that the project should either be saved and finalized, OR should be re-run with different parameters. One parameter, the BioProject ID of the project to check on.
* **`runit`**: Initializes the processing pipeline for a single project and starts the pipeline. This is generally used to start the pipeline *for the first time*, because it creates a new directory for the project and pulls in all the necessary pipeline code. It will throw a warning if the project has been initialized before, but you can safely proceed if the previous run has been removed and you really do want to start over completely. One parameter: BioProject ID of a project for which the metadata is already in our database. (This will retrieve and process the actual FASTQ files.)
* **`discard`**: Throws out any computational results from a single project and records in the database that the project should not be re-attempted. Will also prompt you for a brief explanation of why it should be skipped. The only parameter is the BioProject ID of the project to be thrown out.
* **`again`**: Submits a new slurm job to restart the snakemake pipeline for a single project. *This command assumes the pipeline has already been configured.* Used mostly for situations in which a project stalled for reasons that have been remediated. One parameter, the BioProject ID of the project to be restarted.

### Data ingest commands
(Likely only run once or twice for each export of data from NCBI.)

* **`xml`**: Loads basic sample data exported from the BioProject portal (see "Downloading metadata from NCBI" below for details) and records it in the database. Two required parameters:
  * **taxon** – the NCBI taxon ID used in the search (e.g. txid408170)
  * **filename** - the path to the XML file to be parsed.
* **`tags`**: Similar to the `xml` command and deals with the same input file, but this records all of the tags attached to each sample, rather than their basic metadata. This will likely take much longer. Two required parameters:
  * **taxon** – the NCBI taxon ID used in the search (e.g. txid408170)
  * **filename** - the path to the XML file to be parsed.
* **`runs`**: Queries the compendium database for samples that have an SRS (sample) number, but not an SRR (run) number. This list is then sent to the NCBI eUtils API to retrieve the runs. The only parameter is a limit on how many samples to evaluate (default 2000).

## Downloading metadata from NCBI
We currently extract relevant samples from search results on [the BioSample website](https://www.ncbi.nlm.nih.gov/biosample) using this query:
```
txid408170[Organism:noexp] AND ("public"[filter] AND "biosample sra"[filter])
```
Once the results are displayed, select "Send to" > "File" with the "Full XML (text)" format. This will probably take a while.

## Project states
The `status` table tracks the state of all projects that have been referred to by the CLI in some way. (Projects that had their metadata loaded, but were never processed, are not there.) The possible conditions are:

* `accession_list_created`: A directory has been created for the project and a list of its samples has been added.
* `running`: A job has been submitted to SLURM to be processed
* `to_re_run`: A paired-end project has been flagged for getting re-run as single-end. (This changes to "running" once it's actually been submitted, which happens almost immediately.)
* `failed`: *Terminal status.* Indicates the project was discarded.
* `complete`: Project is done and had acceptable results that were successfully loaded into the database. **Not** a terminal status.
* `archived`: The project's results were stored in a tarball.
* `done`: *Terminal status.* Indicates the project's results were loaded into the database, its files archived, and its other files cleaned up and deleted.

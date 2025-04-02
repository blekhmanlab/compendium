"""
This module assists in the generation and maintenance of a
database containing microbial ecology data from human microbiome
samples.
"""
import click

from datetime import datetime

import config
import db
import projects
import management

@click.group()
def cli():
    pass

@cli.group()
def compendium():
    pass

@cli.group()
def project():
    pass

@compendium.command()
@click.option('--todo', default=2000, help='Number of samples to annotate in this run', show_default=True)
@click.option('--perquery', default=80,
    help='Number of samples to request in each web request. Mostly limited by URL length.',
    show_default=True)
def runs(todo, perquery):
    """Queries the compendium database for samples that have an SRS (sample) number, but not an
    SRR (run) number. This list is then sent to the NCBI eUtils API to retrieve the runs. This
    is required for downloading the raw data.
    """
    db.find_runs(todo, per_query=perquery)

@compendium.command(short_help='Infer a hypervariable region used by each project')
@click.option('--todo', default=100, help='Number of projects to annotate in this run', show_default=True)
def asvs(todo):
    """Runs a heuristic process for inferring which hypervariable regions were
    targeted in an amplicon sequencing project. Records data in the "projects" table.
    """
    db.find_asv_data(todo)

@compendium.command()
@click.argument('taxid')
@click.argument('file')
@click.option('--tags/--skip-tags', default=True, show_default=True,
    help='Indicates whether sample metadata, included in the XML file as key/value pairs,\
    should be recorded in the database.')
def xml(taxid, file, tags):
    """Parse exported BioSample search results and load sample data into the database.

    TAXID is the NCBI taxon ID associated with your samples (e.g. txid408170)
    FILE is the relative path to the XML file to be loaded (e.g. txid408170.xml)
    """
    db.load_xml(taxid, file, save_samples=True, save_tags=tags)

@project.command()
@click.argument('projectid')
def runit(projectid):
    """Process a single project for the first time.
    Initializes the processing pipeline for a single project and starts the pipeline.
    This is generally used to start the pipeline *for the first time*, because it
    creates a new directory for the project and pulls in all the necessary pipeline
    code. It will throw a warning if the project has been initialized before, but
    you can safely proceed if the previous run has been removed and you really do
    want to start over completely.

    This will retrieve and process the actual FASTQ files.

    PROJECTID is a BioProject ID (e.g. PRJNA12345) of a project for which the
        metadata is already in our database.
    """
    proj = projects.Project(project)
    connection = db.Connection()
    proj.initialize_pipeline(connection)
    proj.RUN(connection)

@project.command()
@click.argument('projectid')
def discard(projectid):
    """Throws out any computational results from a single project and records in the
    database that the project should not be re-attempted. Will also prompt you for
    a brief explanation of why it should be skipped. The only parameter is the
    BioProject ID of the project to be thrown out.

    PROJECTID is a BioProject ID (e.g. PRJNA12345)
    """
    proj = projects.Project(projectid)

    confirm = input(f'Really discard project {projectid}? (y/n) ')
    if confirm != 'y':
        click.secho('User input was not "y"; skipping.', fg='red')
        exit(0)

    REASON = input('Provide reason for DB: ')

    proj.errors.append(REASON)
    connection = db.Connection()
    proj.Discard(connection)

@project.command()
@click.argument('projectid')
def again(projectid):
    """Submits a new slurm job to restart the snakemake pipeline for a single
    project. *This command assumes the pipeline has already been configured.*
    Used mostly for situations in which a project stalled for reasons that
    have been remediated. One parameter, the BioProject ID of the project to
    be restarted.

    PROJECTID is a BioProject ID (e.g. PRJNA12345)
    """
    proj = projects.Project(projectid)
    connection = db.Connection()
    proj.RUN(connection)

@project.command()
@click.argument('projectid')
def status(projectid):
    """Retrieves the pipeline progress of a single project and prints a
    report for the user.

    PROJECTID is a BioProject ID (e.g. PRJNA12345)
    """
    proj = projects.Project(projectid)
    if proj.check_if_done(): # true if it's complete
        proj.Load_results_summary()
        proj.print_errors()
    else:
        proj.Report_progress()

@project.command()
@click.argument('projectid')
def eval(projectid):
    """Evaluate the results of a single study. If it's completed the pipeline,
    it will evaluate the results and prompt the user to confirm that the
    project should either be saved and finalized, OR should be re-run with
    different parameters. One parameter, the BioProject ID of the project to
    check on.

    PROJECTID is a BioProject ID (e.g. PRJNA12345)
    """
    proj = projects.Project(projectid)
    if proj.Report_progress(): # true if it's complete
        proj.Load_results_summary()
        proj.print_errors()

        connection = db.Connection()
        proj.REACT(connection)

@compendium.command()
def summary():
    """Summarize the content of the compendium.
    """
    connection = db.Connection()
    management.print_compendium_summary(connection)

@compendium.command()
def status():
    """Summarize the status of any projects with steps
    remaining in their processing pipeline.
    """
    connection = db.Connection()
    current = management.determine_projects(connection)
    management.print_projects_summary(*current)

@compendium.command()
def FORWARD():
    """Interactive process for evaluating all currently pending
    projects. Prompts the user to decide how to deal with
    results.
    """
    connection = db.Connection()
    current = management.determine_projects(connection)
    management.print_projects_summary(*current)
    management.advance_projects(*current, connection)

@compendium.command()
def autoforward():
    """Similar to FORWARD, but automatically approves actions that need to be
    taken. If projects are completed, the application will then search for
    new projects to start.

    Unlike the FORWARD command, this launches new projects as
    others are completed.
    """
    connection = db.Connection()
    # Process the existing projects:
    current = management.determine_projects(connection)
    management.print_projects_summary(*current)
    management.advance_projects(*current, connection, auto=True)

    # Trigger new jobs automatically
    done, running, not_done = current # just unpacking
    TOSTART = config.max_projects-len(running+not_done)

    todo = []
    if TOSTART > 0:
        todo = management.find_todo(connection, needed=TOSTART, max_samples=1000)

    now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    click.secho(
        f"{now}: {len(running+not_done)} projects running. Starting {len(todo)} additional projects: {todo}"
    )
    for pid in todo:
        click.secho(f'Launching {pid}', fg='green')
        proj = projects.Project(pid)
        proj.initialize_pipeline(connection)
        proj.RUN(connection)

if __name__ == "__main__":
    cli()

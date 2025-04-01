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

@cli.command()
@click.option('--todo', default=2000, help='Number of samples to annotate in this run')
@click.option('--perquery', default=80, help='Number of samples to request in each web request. Mostly limited by URL length.')
def runs():
    """Sends requests to the NCBI servers to annotate BioSamples with their SRA run accession codes.
    """
    db.find_runs(todo, per_query=perquery)

@cli.command()
@click.option('--todo', default=100, help='Number of projects to annotate in this run')
def asvs():
    """Runs a heuristic process for inferring which hypervariable regions were
    targeted in an amplicon sequencing project.
    """
    db.find_asv_data(100)

@cli.command()
@click.argument('taxid')
@click.argument('file')
def xml():
    """Parse exported BioSample search results. This loads sample data, but
    skips key/value pairs associated with each sample.

    TAXID is the NCBI taxon ID associated with your samples (e.g. txid408170)
    FILE is the relative path to the XML file to be loaded (e.g. txid408170.xml)
    """
    db.load_xml(taxid, file, save_samples=True, save_tags=False)

@cli.command()
@click.argument('taxid')
@click.argument('file')
def tags():
    """Parse exported BioSample search results. This loads key/value pairs associated
    with each sample, but does not populate the "samples" table itself.

    TAXID is the NCBI taxon ID associated with your samples (e.g. txid408170)
    FILE is the relative path to the XML file to be loaded (e.g. txid408170.xml)
    """
    db.load_xml(taxid, file, save_samples=False, save_tags=True)

@cli.command()
@click.argument('projectid')
def runit():
    """Process a single project for the first time.

    PROJECTID is a BioProject ID (e.g. PRJNA12345)
    """
    proj = projects.Project(project)
    connection = db.Connection()
    proj.initialize_pipeline(connection)
    proj.RUN(connection)

@cli.command()
@click.argument('projectid')
def discard():
    """Clean up a failed project.

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

@cli.command()
@click.argument('projectid')
def again():
    """Retry a failed project.

    PROJECTID is a BioProject ID (e.g. PRJNA12345)
    """
    proj = projects.Project(projectid)
    connection = db.Connection()
    proj.RUN(connection)

@cli.command()
@click.argument('projectid')
def status():
    """Check the status of a single project.

    PROJECTID is a BioProject ID (e.g. PRJNA12345)
    """
    proj = projects.Project(projectid)
    if proj.check_if_done(): # true if it's complete
        proj.Load_results_summary()
        proj.print_errors()
    else:
        proj.Report_progress()

@cli.command()
@click.argument('projectid')
def eval():
    """Evaluate the progress of a single project and process
    its results into the database if appropriate.

    PROJECTID is a BioProject ID (e.g. PRJNA12345)
    """
    proj = projects.Project(projectid)
    if proj.Report_progress(): # true if it's complete
        proj.Load_results_summary()
        proj.print_errors()

        connection = db.Connection()
        proj.REACT(connection)

@cli.command()
def compendium():
    """Summarize the content of the compendium.
    """
    connection = db.Connection()
    management.print_compendium_summary(connection)

@cli.command()
def summary():
    """Summarize the status of any projects with steps
    remaining in their processing pipeline.
    """
    connection = db.Connection()
    current = management.determine_projects(connection)
    management.print_projects_summary(*current)

@cli.command()
def FORWARD():
    """Interactive process for evaluating all currently pending
    projects. Prompts the user to decide how to deal with
    results.
    """
    connection = db.Connection()
    current = management.determine_projects(connection)
    management.print_projects_summary(*current)
    management.advance_projects(*current, connection)

@cli.command()
def autoforward():
    """Interactive process for evaluating all currently pending
    projects. Prompts the user to decide how to deal with
    results.

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

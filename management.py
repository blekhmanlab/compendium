"""
This module deals with automation of project-level operations that
didn't fit with the "projects" module. Selecting which projects to process,
for example.
"""
import projects

import click

def determine_projects(connection):
    """
    Fetches a list of projects without a terminal status and determines
    their current state. Returns 3 arrays of projects:
    1) done
    2) running
    3) not_done (neither done nor running)
    """

    todo = connection.read("""
        SELECT project FROM status
        WHERE status NOT IN ('done','failed')
    """)
    if todo is None:
        click.secho('No projects to evaluate. Exiting.', fg='yellow')
        exit(0)

    todo = [x[0] for x in todo]

    done = []
    running = []
    not_done = [] # jobs that aren't done AND aren't running
    for pid in todo:
        proj = projects.Project(pid)
        if proj.check_if_done():
            done.append(proj)
        else:
            if proj.check_if_running():
                running.append(proj)
            else:
                not_done.append(proj)
    return(done, running, not_done)

def print_projects_summary(done, running, not_done):
    """
    Helper function that prints out lists of BioProject IDs according
    to status.
    """
    click.secho('\n===DONE:', fg='green')
    [click.secho(f'   {x}', bg='green') for x in done]

    click.secho('\n===RUNNING:', fg='yellow')
    [click.secho(f'   {x}', bg='yellow') for x in running]

    click.secho('\n===INCOMPLETE:', fg='red')
    [click.secho(f'   {x}', bg='red') for x in not_done]

def advance_projects(done, running, not_done, connection, auto=False):
    """
    Iterates through projects with things that need to be addressed
    and prompts the user to approve the actions. Projects that are still
    running, or that failed for unknown reasons, simply have their status
    printed.
    """
    for proj in done:
        proj.Load_results_summary()
        proj.print_errors()
        proj.REACT(connection)

    # If this is part of an automated process, don't bother printing progress reports:
    if auto:
        return()

    if len(running) > 0:
        click.secho("\n------------\nSome projects are still running:", fg='yellow')
    for proj in running:
        confirm = input('Print next project? ')
        if confirm != 'y':
            click.secho('Response was not "y"; bailing.', fg='red')
            exit(0)
        proj.Report_progress()

    if len(not_done) > 0:
        click.secho("\n------------\nSome projects are incomplete:", fg='yellow')
    for proj in not_done:
        confirm = input('Print next project? ')
        if confirm != 'y':
            click.secho('Response was not "y"; bailing.', fg='red')
            exit(0)
        proj.Report_progress()

def find_todo(connection, needed=1, min_samples=50, max_samples=10000):
    """
    Reviews the database for projects that have not yet been processed and returns
    a single ID
    """
    done_list = connection.read("""
        SELECT project FROM status
    """)
    if done_list is None:
        done = []
    else:
        done = [x[0] for x in done_list]

    click.secho(f'Tracking down {needed} projects')

    todo = connection.read("""
    SELECT project
    FROM (
        SELECT project, COUNT(srr) AS samples
        FROM SAMPLES s
        WHERE srr IS NOT NULL
            AND library_source IN ('GENOMIC','METAGENOMIC')
            AND library_strategy='AMPLICON'
        GROUP BY 1
        ORDER BY 2 DESC
    ) AS samplecounts
    WHERE samples >= ?
        AND samples <= ?
    ORDER BY RANDOM()
    LIMIT ?
    """, (min_samples, max_samples, needed))

    if todo is None:
        click.secho('Did not find any projects to process.', fg='yellow')
        return []
    return [x[0] for x in todo if x not in done]

def print_compendium_summary(connection):
    """
    Reviews the database for projects that have not yet been processed and returns
    a single ID
    """

    counts = connection.read("""
        SELECT COUNT(DISTINCT project), COUNT(srs) FROM samples
    """)
    if counts is None:
        click.secho('No samples found in samples table.', fg='yellow')
        return()
    click.secho(f'Samples table contains:\n{counts[0][1]} samples from\n{counts[0][0]} projects.\n', fg='green')


    counts = connection.read("""
        SELECT COUNT(DISTINCT project), COUNT(sample)
        FROM (
            SELECT s.project, ac.sample, COUNT(ac.entryid)
            FROM asv_counts ac
            LEFT JOIN samples s
                ON ac.sample=s.srr
            GROUP BY 1,2
        )
    """)

    if counts is None:
        click.secho('No projects found in asv_counts table.', fg='yellow')
        return()
    click.secho(f'Results table contains:\n{counts[0][1]} samples from\n{counts[0][0]} projects.\n', fg='green')

    counts = connection.read("""
        SELECT status, COUNT(DISTINCT project)
        FROM status
        GROUP BY 1
        ORDER BY 2 DESC
    """)
    click.secho('PROJECT STATUS FREQUENCY:')
    for entry in counts:
        click.secho(f'{entry[1]} projects: {entry[0]}')

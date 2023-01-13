import projects

def Determine_projects(connection):
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
        print('No projects to evaluate. Exiting.')
        exit(0)

    todo = [x[0] for x in todo]

    done = []
    running = []
    not_done = [] # jobs that aren't done AND aren't running
    for pid in todo:
        proj = projects.Project(pid)
        if proj.Check_if_done():
            done.append(proj)
        else:
            if proj.Check_if_running():
                running.append(proj)
            else:
                not_done.append(proj)
    return(done, running, not_done)

def Print_projects_summary(done, running, not_done):
    print('\n===DONE:')
    [print(f'   {x}') for x in done]

    print('\n===RUNNING:')
    [print(f'   {x}') for x in running]

    print('\n===INCOMPLETE:')
    [print(f'   {x}') for x in not_done]
    print('\n===\n===\n')

def Advance_projects(done, running, not_done, connection, auto=False):
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
        print("\n------------\nSome projects are still running:")
    for proj in running:
        confirm = input('Print next project? ')
        if confirm != 'y':
            print('Response was not "y"; bailing.')
            exit(0)
        proj.Report_progress()

    if len(not_done) > 0:
        print("\n------------\nSome projects are incomplete:")
    for proj in not_done:
        confirm = input('Print next project? ')
        if confirm != 'y':
            print('Response was not "y"; bailing.')
            exit(0)
        proj.Report_progress()

def Find_todo(connection, needed=1, min_samples=50, max_samples=10000):
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

    print(f'Tracking down {needed} projects')

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
        print('Did not find any projects to process!')
        return([])
    return([x[0] for x in todo if x not in done])
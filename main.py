"""
This module assists in the generation and maintenance of a
database containing microbial ecology data from human microbiome
samples.
"""

from datetime import datetime
import sys # for the command-line params

import config
import db.loader
import db.connector
import projects
import management

if __name__ == "__main__":
    # only command-line param is how many to do in this session
    if sys.argv[1] == 'runs':
        TODO = 2000 if len(sys.argv) < 3 else sys.argv[2]
        db.loader.find_runs(TODO, per_query=80)
    elif sys.argv[1] == 'xml':
        if len(sys.argv) < 4:
            print('The "xml" command requires two parameters: a taxon ID (e.g. txid408170) and the name of the file.')
            exit(1)
        db.loader.load_xml(sys.argv[2], sys.argv[3],
            save_samples=False, save_tags=True)
    elif sys.argv[1] == 'runit':
        if len(sys.argv) < 3:
            exit(1)
        proj = projects.Project(sys.argv[2])
        connection = db.connector.Connection()
        proj.Initialize_pipeline(connection)
        proj.RUN(connection)
    elif sys.argv[1] == 'discard':
        if len(sys.argv) < 3:
            exit(1)
        proj = projects.Project(sys.argv[2])

        confirm = input(f'Really discard project {sys.argv[2]}? (y/n) ')
        if confirm != 'y':
            print('User input was not "y"; skipping.')
            exit(0)
        if len(sys.argv) < 4:
            REASON = input('Provide reason for DB: ')
        else:
            REASON = sys.argv[3]
        proj.errors.append(REASON)
        connection = db.connector.Connection()
        proj.Discard(connection)
    elif sys.argv[1] == 'again':
        if len(sys.argv) < 3:
            exit(1)
        proj = projects.Project(sys.argv[2])
        connection = db.connector.Connection()
        proj.RUN(connection)
    elif sys.argv[1] == 'status':
        if len(sys.argv) < 3:
            exit(1)
        PID = sys.argv[2]

        proj = projects.Project(PID)
        if proj.Check_if_done(): # true if it's complete
            proj.Load_results_summary()
            proj.print_errors()
        else:
            proj.Report_progress()

    elif sys.argv[1] == 'eval':
        if len(sys.argv) < 3:
            exit(1)
        PID = sys.argv[2]

        proj = projects.Project(PID)
        if not proj.Check_if_done(): # true if it's complete
            proj.Report_progress()
            exit(0)

        proj.Load_results_summary()
        proj.print_errors()

        connection = db.connector.Connection()
        proj.REACT(connection)
    elif sys.argv[1] == 'compendium':
        connection = db.connector.Connection()
        management.Print_compendium_summary(connection)
    elif sys.argv[1] == 'summary':
        connection = db.connector.Connection()
        current = management.Determine_projects(connection)
        management.Print_projects_summary(*current)
    elif sys.argv[1] == 'FORWARD':
        connection = db.connector.Connection()
        current = management.Determine_projects(connection)
        management.Print_projects_summary(*current)
        management.Advance_projects(*current, connection)
    elif sys.argv[1] == 'autoforward':
        connection = db.connector.Connection()
        # Process the existing projects:
        current = management.Determine_projects(connection)
        management.Print_projects_summary(*current)
        management.Advance_projects(*current, connection, auto=True)

        # Trigger new jobs automatically
        done, running, not_done = current # just unpacking
        to_start = config.max_projects-len(running+not_done)

        todo = []
        if to_start > 0:
            todo = management.Find_todo(connection, needed=to_start, max_samples=1000)

        now = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
        print(f'{now}: {len(running+not_done)} projects running. Starting {len(todo)} additional projects: {todo}')
        for pid in todo:
            print(f'Launching {pid}')
            proj = projects.Project(pid)
            proj.Initialize_pipeline(connection)
            proj.RUN(connection)

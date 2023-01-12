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
        todo = 2000 if len(sys.argv) < 3 else sys.argv[2]
        db.loader.find_runs(todo, per_query=80)
    elif sys.argv[1] == 'results':
        projects.Process_summary('PRJNA842201')
    elif sys.argv[1] == 'tags':
        db.loader.load_xml('txid408170', 'txid408170.221220.xml',
            save_samples=False, save_tags=True)
    elif sys.argv[1] == 'lists':
        minsamples = 50 if len(sys.argv) < 3 else int(sys.argv[2])
        maxsamples = 50000 if len(sys.argv) < 4 else int(sys.argv[3])
        db.loader.write_lists(minsamples, maxsamples)
    elif sys.argv[1] == 'runit':
        if len(sys.argv) < 3:
            exit(1)
        proj = projects.Project(sys.argv[2])
        connection = db.connector.Connection()
        proj.Initialize_pipeline(connection)
        proj.RUN(connection)
    elif sys.argv[1] == 'again':
        if len(sys.argv) < 3:
            exit(1)
        proj = projects.Project(sys.argv[2])
        connection = db.connector.Connection()
        proj.RUN(connection)
    elif sys.argv[1] == 'status':
        if len(sys.argv) < 3:
            exit(1)
        pid = sys.argv[2]

        proj = projects.Project(pid)
        if proj.Check_if_done(): # true if it's complete
            proj.Load_results_summary()
            proj.print_errors()
        else:
            proj.Report_progress()

    elif sys.argv[1] == 'eval':
        if len(sys.argv) < 3:
            exit(1)
        pid = sys.argv[2]

        proj = projects.Project(pid)
        if not proj.Check_if_done(): # true if it's complete
            proj.Report_progress()
            exit(0)

        proj.Load_results_summary()
        proj.print_errors()

        connection = db.connector.Connection()
        proj.REACT(connection)
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

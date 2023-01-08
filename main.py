import sys # for the command-line params

import db.loader
import db.connector
import results.parsing
import projects.parsing

if __name__ == "__main__":
    # only command-line param is how many to do in this session
    if sys.argv[1] == 'runs':
        todo = 2000 if len(sys.argv) < 3 else sys.argv[2]
        db.loader.find_runs(todo, per_query=80)
    elif sys.argv[1] == 'results':
        projects.parsing.Process_summary('PRJNA842201')
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
        proj = projects.parsing.Project(sys.argv[2])
        connection = db.connector.Connection()
        proj.Initialize_pipeline(connection)
        proj.RUN()
    elif sys.argv[1] == 'again':
        if len(sys.argv) < 3:
            exit(1)
        proj = projects.parsing.Project(sys.argv[2])
        # skip initialization if it's just to restart snakemake
        proj.RUN()
    elif sys.argv[1] == 'status':
        if len(sys.argv) < 3:
            exit(1)
        pid = sys.argv[2]

        proj = projects.parsing.Project(pid)
        if proj.Check_if_done(): # true if it's complete
            proj.Load_results_summary()
            proj.print_errors()

    elif sys.argv[1] == 'eval':
        if len(sys.argv) < 3:
            exit(1)
        pid = sys.argv[2]

        proj = projects.parsing.Project(pid)
        if not proj.Check_if_done(): # true if it's complete
            proj.Report_progress()
            exit(0)

        proj.Load_results_summary()
        proj.print_errors()

        connection = db.connector.Connection()
        proj.REACT(connection)
    elif sys.argv[1] == 'FORWARD':
        connection = db.connector.Connection()
        todo = connection.read("""
            SELECT project FROM status
            WHERE status NOT IN ('done','failed')
        )
        """
        todo = [x[0] for x in todo]
        print(todo)

        not_done = []
        for pid in todo:
            proj = projects.parsing.Project(pid)
            if not proj.Check_if_done():
                not_done.append(pid)
                continue
            proj.Load_results_summary()
            proj.print_errors()
            connection = db.connector.Connection()
            proj.REACT(connection)
        if len(not_done) > 0:
            print("\n------------\nSome projects are incomplete:")
        for pid not_done:
            confirm = input('Print next project?')
            if confirm != 'y':
                print('Response was not "y"; bailing.')
                exit(0)
            proj = projects.parsing.Project(pid)
            proj.Report_progress()

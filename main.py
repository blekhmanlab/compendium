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
        # id = 'PRJNA842201'
        if len(sys.argv) < 3:
            print("NO PROJECT ID PASSED TO COMMAND. Exiting.")
            exit(1)
        proj = projects.parsing.Project(sys.argv[2])
        connection = db.connector.Connection()
        proj.Initialize_pipeline(connection)
        proj.RUN()
    elif sys.argv[1] == 'again':
        if len(sys.argv) < 3:
            print("NO PROJECT ID PASSED TO COMMAND. Exiting.")
            exit(1)
        proj = projects.parsing.Project(sys.argv[2])
        # skip initialization if it's just to restart snakemake
        proj.RUN()

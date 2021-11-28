# This script processes a CSV maintained on the compute cluster
# called "activity.csv" that contains a record of every time a
# new compute job stops or starts. It detects when a study didn't
# complete the pipeline and requires manual intervention.
import csv
from collections import defaultdict

tasks = [
    'download',
    'unpack',
    'dada',
    'trimdada',
    'archive'
]
# some projects are going to have problems we don't want to address
# until later, so leave them out of the list:
with open('/home/blekhman/shared/compendium/code/to_ignore.csv') as f:
    ignore = f.read().splitlines()

progress = defaultdict(lambda: {
    'download': [],
    'unpack': [],
    'dada': [],
    'trimdada': [],
    'archive': []
})

with open('/home/blekhman/shared/compendium/activity.csv','r') as f:
    reader = csv.reader(f)
    for line in reader:
        if len(line) != 3:
            print("!!! Unexpected formatting:") # if two writes happen simultaneously?
            print(line)
            print("Exiting.")
            exit(1)
        project, task, event = line
        progress[project][task].append(event)
to_archive = []
for study in progress.keys():
    if study in ignore: continue
    for task in tasks:
        if 'end' not in progress[study][task]:
            if 'start' not in progress[study][task]:
                if task=='archive':
                    to_archive.append(study)
                elif task == 'trimdada':
                    continue # Not every project needs this
                else:
                    print(f'{study} never started {task}')
                break
            else:
                # If a task started but didn't finish
                if task == 'dada':
                    # if the dada task didn't finish, but re-running
                    # it as trimdada DID finish, chill out
                    if 'end' in progress[study]['trimdada']:
                        continue
                print(f'{study} started {task} but did not finish')
                break
if len(to_archive) > 0:
    print('\n--\nTO ARCHIVE:')
    [print(x) for x in to_archive]

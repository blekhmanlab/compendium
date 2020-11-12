# This script processes a CSV maintained on the compute cluster
# called "activity.csv" that contains a record of every time a
# new compute job stops or starts. It detects when a study didn't
# complete the pipeline and requires manual intervention.
import csv
from collections import defaultdict

progress = defaultdict(lambda: {
    'download': [],
    'trim': [],
    'dada': [],
    'archive': []
})

with open('../activity.csv','r') as f:
    reader = csv.reader(f)
    for line in reader:
        if len(line) != 3:
            print("!!! Unexpected formatting:") # if two writes happen simultaneously?
            print(line)
            print("Exiting.")
            exit(1)
        project, task, event = line
        progress[project][task].append(event)
done = []
for study in progress.keys():
    next_study = False # whether to move on
    for task in ['download','trim','dada','archive']:
        if 'end' not in progress[study][task]:
            if 'start' not in progress[study][task]:
                print(f'{study} never started {task}')
                break
            else:
                print(f'{study} started {task} but did not finish')
                break
    if 'end' in progress[study]['archive']:
        done.append(study)

print("----\nDONE:\n---")
print(done)

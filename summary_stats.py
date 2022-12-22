import config

class Project:
    def __init__(self, name, samples):
        self.id = name
        self.samples = samples

class Sample:
    def __init__(self, data):
        self.srr = data.get('srr')
        self.input = int(data.get('dinput'))
        self.filter = int(data.get('filter'))
        self.forward = int(data.get('forwd'))
        self.length_filter = int(data.get('length'))
        self.nonchim = int(data.get('nonchim'))
        self._check_chimera()

        self.is_paired = 'revse' in data.keys()
        if self.is_paired:
            self.reverse = int(data.get('revse'))
            self.merged = int(data.get('merged'))
            self._check_merged()
            
    def _check_chimera(self):
        """Calculates the proportion of chimeric reads found
        in the reads that were not lost in some other way."""
        self.chimera_percent = 1-(self.nonchim / self.length_filter)
        self.chimeric_warning = self.chimera_percent > config.chimera_worrisome
        self.chimeric_error = self.chimera_percent > config.chimera_too_high
    
    def _check_merged(self):
        """Calculates the proportion of forward reads that were merged
        with a reverse read."""
        self.merged_percent = self.merged / self.forward
        self.merged_warning = self.merged_percent < config.merged_worrisome
        self.merged_error = self.merged_percent < config.merged_too_low


def load_summary(project):
    """Loads a tab-delimited file output by the DADA2 script in which each row
    describes the read counts for a single sample at various steps."""

    samples = []

    with open(f'results/{project}/summary.tsv', 'r') as file:
        headers = file.readline().split('\t')[1:] # first entry is blank
        headers[-1] = headers[-1][:-1] # strip newline
        headers = ['srr'] + headers
        for line in file:
            line = line.split('\t')
            line[-1] = line[-1][:-1] # newline
            data = {}
            for i, entry in enumerate(line):
                data[headers[i]] = entry
            samples.append(Sample(data))
    return(samples)

def Process_summary(project):
    proj = Project(project, load_summary(project))
    print(proj.samples)
    
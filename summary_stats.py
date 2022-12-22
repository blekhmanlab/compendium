import config

class Project:
    def __init__(self, name, samples):
        self.id = name
        self.samples = samples
        self.sample_count = len(self.samples)

        self._check_chimera()
        self._check_merged()
        self._evaluate_flags()

    def _check_chimera(self):
        """Determines what proportion of samples have excessive
        chimeric reads"""
        self.chimeric_warn = 0
        self.chimeric_error = 0
        for sample in self.samples:
            if sample.chimeric_warn:
                self.chimeric_warn += 1
            if sample.chimeric_error:
                self.chimeric_error += 1

    def _check_merged(self):
        """Determines what proportion of samples had an unacceptably
        low proportion of reads merged"""
        merged_warn = 0
        merged_error = 0
        self.paired = True

        for sample in self.samples:
            if not sample.is_paired:
                self.paired = False
                self.merged_warn = None
                self.merged_error = None
                break

            if sample.merged_warn:
                self.merged_warn += 1
            if sample.merged_error:
                self.merged_error += 1
        self.merged_warn = merged_warn / len(self.samples)
        self.merged_error = merged_error / len(self.samples)

    def _evaluate_flags(self):
        """Checks project stats against configured thresholds."""
        self.merged_warn_rate = self.merged_warn / len(self.samples)
        self.merged_error_rate = self.merged_error / len(self.samples)

    def __repr__(self):
        return f'(Project {self.id})'
    def __str__(self):
        return f'Project {self.id}'


class Sample:
    def __init__(self, data):
        self.srr = data.get('srr')[:-8]
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
        self.chimeric_warn = self.chimera_percent > config.chimera_worrisome
        self.chimeric_error = self.chimera_percent > config.chimera_error

    def _check_merged(self):
        """Calculates the proportion of forward reads that were merged
        with a reverse read."""
        self.merged_percent = self.merged / self.forward
        self.merged_warn = self.merged_percent < config.merged_worrisome
        self.merged_error = self.merged_percent < config.merged_error

    def __repr__(self):
        return f'(Sample {self.srr})'
    def __str__(self):
        return f'Sample {self.srr}'


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

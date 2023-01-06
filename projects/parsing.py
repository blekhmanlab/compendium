from datetime import datetime
import glob
import os
import shutil

import config

class Project:
    def __init__(self, name):
        """
        Each Project instance stores metadata about a single BioProject and
        """
        self.id = name # BioProject ID, e.g. PRJNA12345
        self.samples = [] # each item is instance of Sample. Generated from parsing results.
        self.paired = None # is the project's data paired-end?
        self.discard = False # should the project be discarded?
        self.re_run = False # should the project be re-run?
        self.errors = [] # list of issues discovered in the pipeline results

    ##########################
    # Methods for initializing the pipeline
    ##########################
    def _generate_accession_file(self, connection):
        """
        Fetches a list of BioProjects from the local database and generates a list
        of samples for each project.

        Inputs:
            - connection: An instance of type db.connector.Connection
        """
        if self.id is None or len(self.id) < 3:
            raise(Exception(f'Project ID "{self.id}" is not valid. Bailing.'))

        samples = connection.read("""
            SELECT s.srr
            FROM SAMPLES s
            WHERE srr IS NOT NULL
                AND library_source IN ('GENOMIC','METAGENOMIC')
                AND library_strategy='AMPLICON'
            AND project=?""", (self.id,))
        samples = [x[0] for x in samples]

        if os.path.exists(f'{self.id}/SraAccList.txt'):
            print(f'Project already recorded: {self.id}')
            return False

        with open(f'{self.id}/SraAccList.txt','w') as f:
            for sample in samples:
                f.write(f'{sample}\n')
        return True

    def Initialize_pipeline(self, connection):
        """
        Runs a project through the pipeline for the first time.
        """
        os.system(f"git clone --single-branch --depth 1 {config.snakemake_git} {self.id}")
        self._generate_accession_file(connection)

    def RUN(self):
        """
        Starts the pipeline!!!
        """
        timestamp = int(round(datetime.now().timestamp()))
        # go to the project's directory, start the pipeline, then return to
        # wherever we were before
        original = os.getcwd()
        os.chdir(self.id)
        os.system(f'sbatch --job-name={self.id[5:]} -o {self.id}.{timestamp}.log run_snakemake.slurm')
        os.chdir(original)


    ##########################
    # Methods for evaluating processing results
    ##########################
    def Load_results_summary(self):
        """Loads a tab-delimited file output by the DADA2 script in which each row
        describes the read counts for a single sample at various steps.

        Side effects:
        - populates the self.samples list with items of class Sample
        - sets the self.sample_count value
        - performs validation on dada2 results to check for problems with read counts
        - Sets the self.discard and self.re_run flags
        """

        with open(f'{self.id}/summary.tsv', 'r') as file:
            headers = file.readline().split('\t')[1:] # first entry is blank
            headers[-1] = headers[-1][:-1] # strip newline
            headers = ['srr'] + headers
            for line in file:
                line = line.split('\t')
                line[-1] = line[-1][:-1] # newline
                data = {}
                for i, entry in enumerate(line):
                    data[headers[i]] = entry
                self.samples.append(Sample(data))

        self.sample_count = len(self.samples)
        self._check_chimera()
        self._check_merged()
        self._evaluate_flags()

    def _check_chimera(self):
        """Determines what proportion of samples have excessive
        chimeric reads"""
        chimeric_warn = 0
        chimeric_error = 0
        for sample in self.samples:
            if sample.chimeric_warn:
                chimeric_warn += 1
            if sample.chimeric_error:
                chimeric_error += 1
        self.chimeric_warn = chimeric_warn / len(self.samples)
        self.chimeric_error = chimeric_error / len(self.samples)

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
        if self.merged_warn > config.project_merged_worrisome:
            self.re_run = True
            self.errors.append(f'{int(self.merged_warn*100)}% of samples had warning for merged read count.')
        if self.merged_error > config.project_merged_error:
            self.re_run = True
            self.errors.append(f'{int(self.merged_error*100)}% of samples had ERROR for merged read count.')

        if self.chimeric_warn > config.project_chimera_worrisome:
            self.discard = True
            self.errors.append(f'{int(self.chimeric_warn*100)}% of samples had warning for chimeric read count.')
        if self.chimeric_error > config.project_chimera_error:
            self.discard = True
            self.errors.append(f'{int(self.chimeric_error*100)}% of samples had ERROR for chimeric read count.')

    # NOTE: THIS METHOD DELETES FILES
    def Rerun_as_single_end(self):
        """When a paired-end dataset should be re-evaluated without the reverse reads."""
        if not self.paired:
            raise(Exception('Cannot re-run project as single-end; it wasnt paired-end to begin with.'))

        self._remove_previous_dada()
        self._remove_reverse_reads()
        #TODO: actually re-run pipeline

    # NOTE: THIS METHOD DELETES FILES
    def _remove_reverse_reads(self):
        """Deletes all reverse read files extracted from a project's SRA file"""

        files = glob.glob(f'{self.id}/fastq/*_2.fastq')
        for f in files:
            try:
                os.remove(f)
            except OSError as e:
                print(f'Error deleting {f}: {e.strerror}')

    # NOTE: THIS METHOD DELETES FILES
    def _remove_previous_dada(self):
        """Deletes files created by the R scripts of a previous
        processing attempt for a single project. This does NOT
        remove files from the fastq/ directory, which are not
        created by DADA2."""
        if self.id is None or len(self.id) == 0:
            raise(Exception(f'Project ID value is unexpected: "{self.id}"'))

        try:
            shutil.rmtree(f'{self.id}/intermediate')
        except FileNotFoundError:
                pass # no guarantee it was even made
        except OSError as e:
            print(f'Error deleting dir {self.id}/intermediate: {e.strerror}')

        files = [
            'filtered_out.rds',
            'forward_error_model.pdf',
            'reverse_error_model.pdf',
            'err_forward_reads.rds',
            'err_reverse_reads.rds',
            'ASV.tsv', 'asv.rds',
            'ASVs.fa','ASVs_counts.tsv',
            'ASVs_taxonomy.tsv'
        ]
        for f in files:
            try:
                os.remove(f)
            except FileNotFoundError:
                pass # if it's gone, it's fine
            except OSError as e:
                print(f'Error deleting {f}: {e.strerror}')
        # Don't delete the old summary file, just scoot it elsewhere
        if os.path.exists(f'{self.id}/previous_summary.tsv'):
            os.rename(f'{self.id}/previous_summary.tsv', f'{self.id}/previous_previous_summary.tsv')
        if os.path.exists(f'{self.id}/summary.tsv'):
            os.rename(f'{self.id}/summary.tsv', f'{self.id}/previous_summary.tsv')

    def print_errors(self):
        if self.paired is None:
            strategy = 'unknown sequencing strategy'
        else:
            strategy = "paired" if self.paired else "single-end"

        print(f'\nPROJECT {self.id} ({strategy})')
        for error in self.errors:
            print(f'  {error}')
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

def Process_summary(project):
    proj = Project(project)
    #proj.Load_results_summary()
    proj.print_errors()
    return(proj)

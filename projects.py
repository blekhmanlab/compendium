from datetime import datetime
import glob
import os
import shutil
import tarfile

import sqlite3

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
        Fetches a BioProject from the local database and generates a list
        of samples in it. Fed into the SRA Toolkit.

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

        if len(samples) == 0:
            raise(Exception(f'Project {self.id} has no samples in the database!'))

        if os.path.exists(f'{self.id}/SraAccList.txt'):
            raise(Exception(f'Project already initialized: {self.id}'))

        with open(f'{self.id}/SraAccList.txt','w') as f:
            for sample in samples:
                f.write(f'{sample}\n')
        return True

    def _set_status(self, connection, status, note1=None, note2=None):
        """Updates the project's information in the table that tracks project progress"""
        connection.write("""
            UPDATE status
            SET status=?
            WHERE project=?
        """, (status, self.id))

        if note1 is not None:
            connection.write("""
                UPDATE status
                SET note1=?
                WHERE project=?
            """, (note1, self.id))
        if note2 is not None:
            connection.write("""
                UPDATE status
                SET note2=?
                WHERE project=?
            """, (note2, self.id))

    def Initialize_pipeline(self, connection):
        """
        Runs a project through the pipeline for the first time.
        """
        try:
            connection.write('INSERT INTO status (project, status) VALUES (?, ?);', (self.id, 'initialized'))
        except sqlite3.IntegrityError as e:
            print(f'Encountered error writing new project entry to database, likely because project already exists.')
            confirm = input('Continue? (y/n) ')
            if confirm != 'y':
                print('Input was not "y"; bailing.')
                exit(0)
        x = os.system(f"git clone --single-branch --depth 1 {config.snakemake_git} {self.id}")
        if x != 0:
            raise(Exception(f'Call to git returned non-zero exit code {x}'))
        self._generate_accession_file(connection)
        self._set_status(connection, 'accession_list_created')

    def RUN(self):
        """
        Starts the pipeline!!!
        """
        timestamp = int(round(datetime.now().timestamp()))
        x = os.system(f'sbatch --job-name={self.id} -o {self.id}.{timestamp}.log --chdir={self.id} run_snakemake.slurm')
        if x != 0:
            raise(Exception(f'Call to git returned non-zero exit code {x}'))

    def Check_if_done(self):
        to_check = [
            f'{self.id}/ASVs.fa',
            f'{self.id}/ASVs_counts.tsv',
            f'{self.id}/ASVs_taxonomy.tsv'
        ]
        # return True if all the required files are present, otherwise False
        return(False not in [os.path.exists(x) for x in to_check])

    def Check_if_running(self):
        return(os.path.exists(f'{self.id}/running.txt'))

    def Report_progress(self):
        print(self)
        if self.Check_if_done():
            print('DONE!')
            return(True)

        if self.Check_if_running():
            print('\n===============\nCURRENTLY RUNNING\n===============\n')

        tests = [
            ('Initialization',
            {
                'Directory created': os.path.isdir(self.id),
                'Repository cloned': os.path.isdir(f'{self.id}/workflow'),
                'Accession list created': os.path.exists(f'{self.id}/SraAccList.txt'),
                'Virtual environment created': os.path.isdir(f'{self.id}/venv')
            }),
            ('Pipeline',
            {
                '1/6 Prefetch job started': os.path.exists(f'{self.id}/.snakemake/slurm_logs/rule_sra_prefetch'),
                '2/6 SRA data extraction job started': os.path.exists(f'{self.id}/.snakemake/slurm_logs/rule_sra_to_fastq'),
                '3/6 FASTQ filtering job started': os.path.exists(f'{self.id}/.snakemake/slurm_logs/rule_filter'),
                '4/6 Error modeling job started': os.path.exists(f'{self.id}/.snakemake/slurm_logs/rule_errormodel'),
                '5/6 ASV calculation job started': os.path.exists(f'{self.id}/.snakemake/slurm_logs/rule_make_asv_table'),
                '6/6 Taxonomic assignment job started': os.path.exists(f'{self.id}/.snakemake/slurm_logs/rule_assign_taxonomy')
            }),
            ('Results',
            {
                f'Result file: ASVs.fa': os.path.exists(f'{self.id}/ASVs.fa'),
                f'Result file: ASVs_counts.tsv': os.path.exists(f'{self.id}/ASVs_counts.tsv'),
                f'Result file: ASVs_taxonomy.tsv': os.path.exists(f'{self.id}/ASVs_taxonomy.tsv'),
            })
        ]

        arrow = True # point at the earliest test that fails
        for category in tests:
            print(f'\n======{category[0]}======')
            for string, test in category[1].items():
                print(f"{'✓' if test else 'X'}   {string} {'  <<< XXXXXXX <<<' if arrow and not test else ''}")
                if arrow and not test:
                    arrow = False # only print one arrow

        return(False)

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

        self._check_retained()
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
        warncount = 0
        errorcount = 0
        self.paired = True

        for sample in self.samples:
            if not sample.is_paired:
                self.paired = False
                self.merged_warn = None
                self.merged_error = None
                break

            if sample.merged_warn:
                warncount += 1
            if sample.merged_error:
                errorcount += 1
        self.merged_warn = warncount / len(self.samples)
        self.merged_error = errorcount / len(self.samples)

    def _check_retained(self):
        """Determines what proportion of samples had an unacceptably
        low proportion of reads merged"""
        warncount = 0
        errorcount = 0

        for sample in self.samples:
            if sample.retained_warn:
                warncount += 1
            if sample.retained_error:
                errorcount += 1
        self.retained_warn = warncount / len(self.samples)
        self.retained_error = errorcount / len(self.samples)

    def _evaluate_flags(self):
        """Checks project stats against configured thresholds."""
        if self.merged_warn > config.project_merged_worrisome:
            self.re_run = True
            self.errors.append(f'{int(self.merged_warn*100)}% of samples had warning for merged read count.')
        if self.merged_error > config.project_merged_error:
            self.re_run = True
            self.errors.append(f'{int(self.merged_error*100)}% of samples had ERROR for merged read count.')

        # Don't bother checking the percentage of reads retained if the project
        # is going to be re-run as single-ended anyway. We could probably catch
        # a few doomed projects earlier if we check it here, but the logic is more
        # complicated than we should bother with right now.
        if self.re_run:
            return

        # Read retention
        if self.retained_warn > config.project_retained_worrisome:
            self.discard = True
            self.errors.append(f'{int(self.retained_warn*100)}% of samples had warning for reads retained.')
        if self.retained_error > config.project_retained_error:
            self.discard = True
            self.errors.append(f'{int(self.retained_error*100)}% of samples had ERROR for reads retained.')

        # Chimeric reads
        if self.chimeric_warn > config.project_chimera_worrisome:
            self.discard = True
            self.errors.append(f'{int(self.chimeric_warn*100)}% of samples had warning for chimeric read count.')
        if self.chimeric_error > config.project_chimera_error:
            self.discard = True
            self.errors.append(f'{int(self.chimeric_error*100)}% of samples had ERROR for chimeric read count.')

    # NOTE: THIS METHOD DELETES FILES AND STARTS PIPELINES
    def Rerun_as_single_end(self, connection):
        """When a paired-end dataset should be re-evaluated without the reverse reads."""
        if not self.paired:
            raise(Exception('Cannot re-run project as single-end; it wasnt paired-end to begin with.'))

        self._remove_previous_dada()
        self._remove_reverse_reads()
        connection.write("""
            UPDATE status
            SET rerun_as_single_end=1
            WHERE project=?
        """, (self.id,))

        self.RUN()

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

    def _record_if_paired(self, connection):
        """
        Records whether a project has paired-end data. Note that this overwrites
        previous values set in this column, so if a project is marked as single-end,
        that means the final results were single-end, not that the data was. If
        paired is true AND rerun_as_single_end is true, this means "paired-end data"
        also.
        """

        if self.paired is not None:
            print(f"(Recording that project was {'' if self.paired else 'not '}paired-end.)")
            connection.write("""
                UPDATE status
                SET paired=?
                WHERE project=?
            """, (1 if self.paired else 0, self.id))

    def Discard(self, connection):
        """Removes a project's files and records its status as failed"""
        self._record_if_paired(connection)

        print(f'DELETING PROJECT {self.id}')
        self._set_status(connection, 'failed', ' / '.join(self.errors))
        shutil.rmtree(f'{self.id}')
        return(True)

    ########### helpers for saving results
    def _load_counts(self):
        """Loads a tab-delimited file in which each column
        is a sample and each row is an ASV, with the cells
        indicating read counts."""

        to_write = []

        with open(f'{self.id}/ASVs_counts.tsv', 'r') as file:
            samples = file.readline().split('\t')[1:] # first entry is blank
            samples[-1] = samples[-1][:-1] # strip newline
            for line in file:
                line = line.split('\t')
                line[-1] = line[-1][:-1]
                asv = [line[0]] * len(samples)
                entries = list(zip(samples, asv, line[1:]))
                # convert counts to ints
                entries = [
                    (project, x[0], x[1], int(x[2]))
                    for x in entries
                ]
                to_write += [x for x in entries if x[2] > 0]
        # example entry: ('PRJNA987', 'SRR123', 'ASV_7', 23)
        return(to_write)

    def _load_asv_data(self):
        """Loads a tab-delimited file in which each row is
        a numbered ASV, associated with its inferred taxonomic
        source."""

        seqs = {}
        # Get exact sequences
        with open(f'{self.id}/ASVs.fa', 'r') as file:
            while file:
                asv = file.readline()
                if asv == '': # empty line with no newline
                    break
                if len(asv) > 2:
                    asv = asv[1:-1] # strip off leading '>' and trailing newline
                seq = file.readline()
                if len(seq) > 1:
                    seq = seq[0:-1] # strip trailing newline
                seqs[asv]=[seq]

        assignments = {}
        # The get taxonomic assignments
        with open(f'{self.id}/ASVs_taxonomy.tsv', 'r') as file:
            file.readline() # skip header
            for line in file:
                line = line.split('\t')
                line[-1] = line[-1][:-1]
                assignments[line[0]] = line[1:]

        # example entry: (
        #       ('ASV_1','Bacteria','Bacteroidota','Bacteroidia','Bacteroidales','Bacteroidaceae','Bacteroides'),
        #       ('PRJNA1234', 'ASV_1','CCTACGGG')
        # )
        return([(tuple([asv]+values), (project, asv, seqs[asv])) for asv, values in assignments.items()])

    def Save_results(self, connection):
        """
        Loads the DADA2 results and saves them to the DB
        """
        self._record_if_paired(connection)

        counts = self._load_counts()
        assignments, seqs = self._load_asv_data()

        connection.write('INSERT INTO asv_counts VALUES (?,?,?,?)', counts)
        asv_ids = connection.write("""
            INSERT INTO asv_sequences (project, asv, seq)
            VALUES (?,?,?)
            RETURNING asv, asv_id
        """, seqs)
        ids = {}
        # figure out which ASV ID goes with which ASV we just recorded:
        for asv, asv_id in asv_ids:
            ids[asv] = asv_id
        # each assignment entry has a project-level ASV id (ASV_1, ASV_2, etc), but
        # we want to swap that out for the unique ID assigned by SQLite when we saved the ASV's sequence:
        to_write = [tuple([ids[entry[0]], 'silva_nr99_v138_train_set']+entry[1:]) for entry in assignments]
        asv_ids = connection.write("""
            INSERT INTO asv_assignments
            VALUES (?,?,?,?,?,?,?,?)
        """, to_write)

        self._set_status(connection, 'complete')

        confirm = input('Results recorded. Archive results? ')
        if confirm != 'y':
            print('User input was not "y"; skipping.')
            return

        if not os.path.exists('archive'):
            os.mkdir('archive')

        with tarfile.open(name=f'archive/{self.id}.tar.gz', mode='w:gz') as archive:
            archive.add(f'{self.id}/.snakemake/log')
            archive.add(f'{self.id}/.snakemake/slurm_logs')
            archive.add(f'{self.id}/ASVs_taxonomy.tsv')
            archive.add(f'{self.id}/ASVs.fa')
            archive.add(f'{self.id}/ASVs_counts.tsv')
            # find the log file
            for f in os.listdir(self.id):
                if f.endswith('.log') or f.endswith('.pdf'):
                    archive.add(f'{self.id}/{f}')

        if not os.path.exists(f'archive/{self.id}.tar.gz'):
            raise(Exception(f'Archive of project {self.id} was not found in archive directory.'))
        self._set_status(connection, 'archived')
        confirm = input('Archive created. Delete files? ')
        if confirm != 'y':
            print('User input was not "y"; skipping.')
            return
        shutil.rmtree(f'{self.id}')
        if not os.path.exists(self.id):
            self._set_status(connection, 'done')
        return

    def REACT(self, connection):
        '''Acts on the results of the pipeline completion'''
        if self.discard:
            confirm = input(f'Delete project {self.id}? (y/n) ')
            if confirm != 'y':
                print('Will only delete if user responds "y". Skipping.')
                return
            self.Discard(connection)
            return(True)

        elif self.re_run:
            confirm = input(f'Re-run project {self.id} as single end? (y/n) ')
            if confirm != 'y':
                print('Will only re-run if user responds "y". Skipping.')
                return
            self.Rerun_as_single_end(connection)
            return(True)
        # if we make it to this point, it's good to go!
        print(f'\nProject {self.id} has passed all checks!')
        confirm = input(f'Save results of project {self.id} ({self.sample_count} samples)? (y/n) ')
        if confirm != 'y':
            print('Will only save if user responds "y". Skipping.')
            return()
        self.Save_results(connection)

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

        self.is_paired = 'revse' in data.keys()
        self.reverse = None
        self.merged = None

        self._check_chimera()
        if self.is_paired:
            self.reverse = int(data.get('revse'))
            self.merged = int(data.get('merged'))
            self._check_merged()
        self._check_stages()

    def _check_chimera(self):
        """Calculates the proportion of chimeric reads found
        in the reads that were not lost in some other way."""
        try:
            # This will fail if there's a sample (a control, for example)
            # that has zero reads that cleared the length filter.
            self.chimera_percent = 1-(self.nonchim / self.length_filter)
            self.chimeric_warn = self.chimera_percent > config.chimera_worrisome
            self.chimeric_error = self.chimera_percent > config.chimera_error
        except:
            self.chimeric_error = False
            self.chimeric_warn = False

    def _check_merged(self):
        """Calculates the proportion of forward reads that were merged
        with a reverse read."""
        try:
            # this will fail if there's a sample with zero forward reads that were kept
            self.merged_percent = self.merged / self.forward
            self.merged_warn = self.merged_percent < config.merged_worrisome
            self.merged_error = self.merged_percent < config.merged_error
        except:
            self.merged_error = False
            self.merged_warn = False

    def _check_stages(self):
        self.retained_percent = self.nonchim / self.input
        # TODO: add a check here to measure how many reads dropped
        # out between every stage
        self.retained_warn = self.retained_percent < config.retained_worrisome
        self.retained_error = self.retained_percent < config.retained_error

    def __repr__(self):
        return f'(Sample {self.srr})'
    def __str__(self):
        return f'Sample {self.srr}'

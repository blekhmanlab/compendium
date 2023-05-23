# What identifiers should we send with API requests sent to NCBI?
Tool='compendium'
Email='email_goes_here'
Key='12345' # API key obtained (for free!) from NCBI

# Where should we store the archived data?
archive_path = '/path/to_your/archive/'

# Where is the SQLite database file?
db_path = '/path_to_your/compendium.db'

# How many projects should we try to have running at one time?
max_projects = 8

# Should a user's input be required to reprocess and discard projects?
confirm_destruct = True

# Where does the pipeline code live?
snakemake_git = 'git@github.com:blekhmanlab/snakemake-compendium.git'

############### Web request settings ###############
# When building up metadata for each project, we need to make multiple calls
# to the NCBI's eUtils API. These settings dictate how this happens.
#
# How long should we wait to get a response from a web request before bailing?
timeout = 25
#
# What base URLs should be used when building up queries?
esearch_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?tool={Tool}&email={Email}&api_key={Key}&db=sra&usehistory=y&term=' # pylint: disable=line-too-long
efetch_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?tool={Tool}&email={Email}&db=sra&query_key=1'
#
##############################################################


############### Sample-level settings ###############
# When evaluating the proportion of reads retained through the entire pipeline,
# what thresholds should be used for the "warning" and "error" levels?
retained_worrisome = 0.69
retained_error = 0.59 # MINIMUM
#
# When evaluating the proportion of reads removed from a sample,
# what thresholds should be used for the "warning" and "error" levels?
chimera_worrisome = 0.10
chimera_error = 0.20 # MAXIMUM
#
# When evaluating the proportion of forward reads removed from a
# sample, what thresholds should be used for the "warning" and "error" levels?
merged_worrisome = 0.80
merged_error = 0.65 # MINIMUM
#
##############################################################


############### Project-level settings ###############
#
# When evaluating the proportion of SAMPLES flagged for issues,
# what thresholds should be used for determining whether a project
# should be reprocessed?
#
# (IMPORTANT: The thresholds here work slightly differently
# than the sample-level ones above: There is no "warning" level for projects.
# The "worrisome" threshold here is *what proportion of samples must be marked as
# worrisome for the project to be discarded*. The "error" threshold is what proportion
# of samples must be marked as errors for the project to be discarded.)
project_retained_worrisome = 0.70
project_retained_error = 0.20
#
project_merged_worrisome = 0.4
project_merged_error = 0.15
#
project_chimera_worrisome = 0.4
project_chimera_error = 0.2
#
##############################################################

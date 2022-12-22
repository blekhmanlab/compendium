Tool='compendium'
Email='email_goes_here'
Key='12345'


############### Sample-level settings ###############
#
# When evaluating the proportion of reads removed from a sample,
# what thresholds should be used for the "warning" and "error" levels?
chimera_worrisome = 0.10
chimera_error = 0.60 # MAXIMUM
#
# When evaluating the proportion of forward reads removed from a
# sample, what thresholds should be used for the "warning" and "error" levels?
merged_worrisome = 0.65
merged_error = 0.50 # MINIMUM
##############################################################


############### Project-level settings ###############
#
# When evaluating the proportion of SAMPLES flagged for excessive
# chimeric reads, what thresholds should be used for determining whether
# a project should be reprocessed?
# (IMPORTANT: The thresholds for chimeras and merged reads work slightly differently
# than the sample-level ones above: There is no "warning" level for projects.
# The "worrisome" threshold here is *what proportion of samples must be marked as
# worrisome for the project to be discarded*. The "error" threshold is what proportion
# of samples must be marked as errors for the project to be discarded.)
project_chimera_worrisome = 0.4
project_chimera_error = 0.2
#
# When evaluating the proportion of samples flagged for low
# rates of merged reads, what thresholds should be used?
project_merged_worrisome = 0.4
project_merged_error = 0.15
##############################################################

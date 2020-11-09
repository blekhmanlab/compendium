# https://astrobiomike.github.io/amplicon/dada2_workflow_ex

library(dada2)
packageVersion("dada2")

setwd('/mnt/trimmed')

log <- function(message) print(paste(date(), message))

## first we're setting a few variables we're going to use
# one with all sample names, by scanning our "samples" file we made earlier
samples <- scan("../samples.txt", what="character")

# one holding the file names of all the forward reads
forward_reads <- paste0(samples, ".R1.fq")
# and one with the reverse
reverse_reads <- paste0(samples, ".R2.fq")
# and variables holding file names for the forward and reverse
# filtered reads we're going to generate below
filtered_forward_reads <- paste0("../intermediate/", samples, "_R1_filtered.fastq.gz")
filtered_reverse_reads <- paste0("../intermediate/", samples, "_R2_filtered.fastq.gz")


#########################
# Quality filtering
#########################

# Filter for quality. We do this even if shi7 is used first because
# it does a few DADA2-specific things (throwing out 'N' bases, for example)
# that become important later.
log('Filtering...')
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
    reverse_reads, filtered_reverse_reads,
    rm.phix=TRUE, multithread=8)

#########################
# Building error models
#########################
log('Building forward error model...')
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
log('Building reverse error model...')
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)

log('Plotting error models...')
pdf('../forward_error_model.pdf')
plotErrors(err_forward_reads, nominalQ=TRUE)
dev.off()

pdf('../reverse_error_model.pdf')
plotErrors(err_reverse_reads, nominalQ=TRUE)
dev.off()

#########################
# Dereplicate identical reads
#########################
log('Dereplicating...')
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples

#########################
# Generate count table
#########################
log('Processing forward reads...')
dada_forward <- dada(derep_forward, err=err_forward_reads, multithread=TRUE)
log('Processing reverse reads...')
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, multithread=TRUE)

log('Merging reads...')
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
    derep_reverse, trimOverhang=TRUE, minOverlap=235, maxOverlap=255)

seqtab <- makeSequenceTable(merged_amplicons)
# check for chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=T)

#########################
# Check reads dropped at each step
#########################
getN <- function(x) sum(getUniques(x))

print('Calculating summary stats...')
# making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          dada_r=sapply(dada_reverse, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

# OUTPUT
log('Writing summary output...')
write.table(summary_tab, "../summary.tsv",
            sep="\t", quote=F, col.names=NA)
log('Writing ESV table...')
write.table(seqtab.nochim, "../ASV.tsv",
            sep="\t", quote=F, col.names=NA)
saveRDS(seqtab.nochim, '../asv.rds')
log('DONE!!!')

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
filtered_forward_reads <- paste0("../intermediate/", samples, ".R1.filtered.fastq.gz")
filtered_reverse_reads <- paste0("../intermediate/", samples, ".R2.filtered.fastq.gz")


#########################
# Quality filtering
#########################

# Filter for quality. We do this even if shi7 is used first because
# it does a few DADA2-specific things (throwing out 'N' bases, for example)
# that become important later.
log('Filtering...')
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads,
                              truncQ=0, rm.phix=TRUE, multithread=4)

# saveRDS(filtered_out, '../temp/filtered_out.rds')
# filtered_out <- readRDS('../temp/filtered_out.rds')

#########################
# Building error models
#########################
log('Building forward error model...')
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)
log('Building reverse error model...')
err_reverse_reads <- learnErrors(filtered_reverse_reads, multithread=TRUE)

log('Plotting error models...')
pdf('../temp/forward_error_model.pdf')
plotErrors(err_forward_reads, nominalQ=TRUE)
dev.off()

pdf('../temp/reverse_error_model.pdf')
plotErrors(err_reverse_reads, nominalQ=TRUE)
dev.off()

# saveRDS(err_forward_reads, '../temp/err_forward_reads.rds')
# saveRDS(err_reverse_reads, '../temp/err_reverse_reads.rds')
# err_forward_reads <- readRDS('../temp/err_forward_reads.rds')
# err_reverse_reads <- readRDS('../temp/err_reverse_reads.rds')


#######################
# Grab a sample to process
######################
if(length(samples) > 15) {
  subsamples <- sample(samples, 15)
} else {
  subsamples <- samples
}

subfiltered_forward_reads <- paste0("../intermediate/", subsamples, ".R1.filtered.fastq.gz")
subfiltered_reverse_reads <- paste0("../intermediate/", subsamples, ".R2.filtered.fastq.gz")

# Dereplicate
subderep_forward <- derepFastq(subfiltered_forward_reads, verbose=TRUE)
names(subderep_forward) <- subsamples
subderep_reverse <- derepFastq(subfiltered_reverse_reads, verbose=TRUE)
names(subderep_reverse) <- subsamples
# count table
subdada_forward <- dada(subderep_forward, err=err_forward_reads, multithread=TRUE)
subdada_reverse <- dada(subderep_reverse, err=err_reverse_reads, multithread=TRUE)

getN <- function(x) sum(getUniques(x))
# count the forward reads we classified
dada_f=sapply(subdada_forward, getN)

minoverlap = 0
for(mintest in c(235, 200, 165, 140, 115, 90, 75, 60, 20)) {
  submerged_amplicons <- mergePairs(subdada_forward, subderep_forward, subdada_reverse,
                                 subderep_reverse, trimOverhang=TRUE, minOverlap=mintest)
  
  # count the merged reads we ended up with
  merged=sapply(submerged_amplicons, getN)
  
  if(sum(merged) / sum(dada_f) > 0.65) {
    minoverlap = mintest
    log(paste("Proceeding with minOverlap of", minoverlap, "based on merge rate of", sum(merged) / sum(dada_f)))
    break
  }
  log(paste("Tested", mintest, "but got a merge rate of", (sum(merged) / sum(dada_f))))
}

if(minoverlap == 0) {
  stop("COULDN'T DETERMINE AN OVERLAP SETTING. BAILING.")
}

#########################
# RETURNING TO PROCESSING THE FULL DATASET
#########################
#Dereplicate identical reads
log('Dereplicating...')
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtered_reverse_reads, verbose=TRUE)
names(derep_reverse) <- samples

# saveRDS(derep_forward, '../temp/derep_forward.rds')
# saveRDS(derep_reverse, '../temp/derep_reverse.rds')
# derep_forward <- readRDS('../temp/derep_reverse.rds')
# derep_reverse <- readRDS('../temp/derep_reverse.rds')

#########################
# Generate count table
#########################
log('Processing forward reads...')
dada_forward <- dada(derep_forward, err=err_forward_reads, multithread=TRUE)
log('Processing reverse reads...')
dada_reverse <- dada(derep_reverse, err=err_reverse_reads, multithread=TRUE)

# saveRDS(dada_forward, '../temp/dada_forward.rds')
# saveRDS(dada_reverse, '../temp/dada_reverse.rds')
# dada_forward <- readRDS('../temp/dada_forward.rds')
# dada_reverse <- readRDS('../temp/dada_reverse.rds')

log('Merging reads...')
merged_amplicons <- mergePairs(dada_forward, derep_forward, dada_reverse,
                               derep_reverse, trimOverhang=TRUE, minOverlap=minoverlap)

seqtab <- makeSequenceTable(merged_amplicons)
# check for chimeras
log('Removing bimeras...')
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
log('Writing ASV table...')
write.table(seqtab.nochim, "../ASV.tsv",
            sep="\t", quote=F, col.names=NA)
saveRDS(seqtab.nochim, '../asv.rds')
log('DONE!!!')

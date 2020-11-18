# https://astrobiomike.github.io/amplicon/dada2_workflow_ex

library(dada2)
packageVersion("dada2")

setwd('/mnt/trimmed')

log <- function(message) print(paste(date(), message))

## first we're setting a few variables we're going to use
# one with all sample names, by scanning our "samples" file we made earlier
samples <- scan("../samples.txt", what="character")

# one holding the file names of all the forward reads
forward_reads <- paste0(samples, ".1.fq")
filtered_forward_reads <- paste0("../intermediate/", samples, ".R1.filtered.fastq.gz")


#########################
# Quality filtering
#########################

# Filter for quality. We do this even if shi7 is used first because
# it does a few DADA2-specific things (throwing out 'N' bases, for example)
# that become important later.
log('Filtering...')
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              truncQ=0, rm.phix=TRUE, multithread=4)

# saveRDS(filtered_out, '../temp/filtered_out.rds')
# filtered_out <- readRDS('../temp/filtered_out.rds')

#########################
# Building error models
#########################
log('Building forward error model...')
err_forward_reads <- learnErrors(filtered_forward_reads, multithread=TRUE)

log('Plotting error models...')
pdf('../temp/forward_error_model.pdf')
plotErrors(err_forward_reads, nominalQ=TRUE)
dev.off()

# saveRDS(err_forward_reads, '../temp/err_forward_reads.rds')
# err_forward_reads <- readRDS('../temp/err_forward_reads.rds')

# Dereplicate identical reads
log('Dereplicating...')
derep_forward <- derepFastq(filtered_forward_reads, verbose=TRUE)
names(derep_forward) <- samples # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow

# saveRDS(derep_forward, '../temp/derep_forward.rds')
# derep_forward <- readRDS('../temp/derep_reverse.rds')

#########################
# Generate count table
#########################
log('Processing forward reads...')
dada_forward <- dada(derep_forward, err=err_forward_reads, multithread=TRUE)

# saveRDS(dada_forward, '../temp/dada_forward.rds')
# dada_forward <- readRDS('../temp/dada_forward.rds')

seqtab <- makeSequenceTable(dada_forward)

# Get rid of really short sequences that can't practically be used
# to assign taxonomy:
seqtab.noshort <- seqtab[,nchar(colnames(seqtab)) > 49]
log(paste('Removed',length(colnames(seqtab) - length(colnames(seqtab.noshort)),'ASVs for being too short.')))

# check for chimeras
log('Removing bimeras...')
seqtab.nochim <- removeBimeraDenovo(seqtab.noshort, verbose=T)

#########################
# Check reads dropped at each step
#########################
getN <- function(x) sum(getUniques(x))

print('Calculating summary stats...')
# making a little table
summary_tab <- data.frame(row.names=samples, dada2_input=filtered_out[,1],
                          filtered=filtered_out[,2], dada_f=sapply(dada_forward, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/filtered_out[,1]*100, 1))

# OUTPUT
log('Writing summary output...')
write.table(summary_tab, "../results/summary.tsv",
            sep="\t", quote=F, col.names=NA)
log('Writing ASV table...')
write.table(seqtab.nochim, "../results/ASV.tsv",
            sep="\t", quote=F, col.names=NA)
saveRDS(seqtab.nochim, '../asv.rds')
log('ASVs recorded.')

log('Assigning taxonomy...')
taxa <- assignTaxonomy(seqtab.nochim, "~/silva_nr_v138_train_set.fa.gz", multithread=8, tryRC=T)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
# extract fasta:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
# tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)


## output
log('Writing output files...')
write(asv_fasta, "../results/ASVs.fa")
write.table(asv_tab, "../results/ASVs_counts.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax, "../results/ASVs_taxonomy.tsv",
            sep="\t", quote=F, col.names=NA)

log('DONE!!!')

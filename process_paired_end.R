# https://astrobiomike.github.io/amplicon/dada2_workflow_ex

library(dada2)
packageVersion("dada2")

setwd('/mnt/fastq')

log <- function(message) print(paste(date(), message))

## first we're setting a few variables we're going to use
# one with all sample names, by scanning our "samples" file we made earlier
samples <- scan("../samples.txt", what="character")

# one holding the file names of all the forward reads
forward_reads <- paste0(samples, "_1.fastq")
# and one with the reverse
reverse_reads <- paste0(samples, "_2.fastq")
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

trimLeft = 0
trimRight = 0
args = commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  if (args[1] == 'trim') {
    trimLeft = 24
    trimRight = 24
  }
}
filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              reverse_reads, filtered_reverse_reads,
                              truncQ=0, rm.phix=TRUE, multithread=4,
                              verbose=TRUE, trimLeft=trimLeft, trimRight=trimRight)

# saveRDS(filtered_out, 'filtered_out.rds')
#filtered_out <- readRDS('filtered_out.rds')

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

# saveRDS(err_forward_reads, 'err_forward_reads.rds')
# saveRDS(err_reverse_reads, 'err_reverse_reads.rds')
#err_forward_reads <- readRDS('err_forward_reads.rds')
#err_reverse_reads <- readRDS('err_reverse_reads.rds')

#########################
# Generate count table
#########################
mergers <- vector("list", length(samples))
names(mergers) <- samples

ddF <- vector("list", length(samples))
names(ddF) <- samples
ddR <- vector("list", length(samples))
names(ddR) <- samples

for(sam in samples) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(paste("../intermediate/", sam, ".R1.filtered.fastq.gz", sep=""))
  ddF[[sam]] <- dada(derepF, err=err_forward_reads, multithread=TRUE)
  derepR <- derepFastq(paste("../intermediate/", sam, ".R2.filtered.fastq.gz", sep=""))
  ddR[[sam]] <- dada(derepR, err=err_reverse_reads, multithread=TRUE)
  merger <- mergePairs(ddF[[sam]], derepF, ddR[[sam]], derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

seqtab <- makeSequenceTable(mergers)

# Get rid of really short sequences that can't practically be used
# to assign taxonomy:
seqtab.noshort <- seqtab[,nchar(colnames(seqtab)) > 49]
diff <- length(colnames(seqtab)) - length(colnames(seqtab.noshort))
log(paste('Removed',diff,'ASVs for being too short.'))

# check for chimeras
log('Removing bimeras...')
seqtab.nochim <- removeBimeraDenovo(seqtab.noshort, verbose=T)

#########################
# Check reads dropped at each step
#########################
getN <- function(x) sum(getUniques(x))

print('Calculating summary stats...')
# making a little table
merged_val <- sapply(mergers, getN)
nochim_val <- rowSums(seqtab.nochim)
chim_removed_val <- round(((merged_val-nochim_val)/filtered_out[,1])*100, 1)

summary_tab <- data.frame(row.names=samples, dinput=filtered_out[,1],
                          filter=filtered_out[,2], forwd=sapply(ddF, getN),
                          revse=sapply(ddR, getN), merged=merged_val,
                          nonchim=nochim_val,
                          chim_perc=chim_removed_val,
                          retained_perc=round((nochim_val*100)/filtered_out[,1], 1))

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
taxa <- assignTaxonomy(seqtab.nochim, "/code/silva_nr99_v138_train_set.fa.gz", multithread=8, tryRC=T)

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

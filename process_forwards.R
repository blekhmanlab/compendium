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
filtered_forward_reads <- paste0("../intermediate/", samples, ".R1.filtered.fastq.gz")


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
  }
}

filtered_out <- filterAndTrim(forward_reads, filtered_forward_reads,
                              truncQ=0, rm.phix=TRUE, multithread=4,
                              verbose=TRUE, trimLeft=trimLeft, trimRight=trimRight)

# Then we limit the list of filtered fastq files to include
# only the ones that actually had reads pass the filter:
filtered_forward_reads <- filtered_forward_reads[file.exists(filtered_forward_reads)]

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

#########################
# Generate count table
#########################
ddF <- vector("list", length(samples))
names(ddF) <- samples

for(sam in samples) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(paste("../intermediate/", sam, ".R1.filtered.fastq.gz", sep=""))
  ddF[[sam]] <- dada(derepF, err=err_forward_reads, multithread=TRUE)
}
rm(derepF)

seqtab <- makeSequenceTable(ddF)

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
forwd_val <- sapply(ddF, getN)
nochim_val <- rowSums(seqtab.nochim)
length_val <- rowSums(seqtab.noshort)
chim_removed_val <- round(((length_val-nochim_val)/forwd_val)*100, 1)

summary_tab <- data.frame(dinput=filtered_out[,1],
                          filter=filtered_out[,2], forwd=forwd_val,
                          length=length_val,
                          nonchim=nochim_val,
                          chim_perc=chim_removed_val,
                          retained_perc=round((100*nochim_val)/filtered_out[,1], 1))

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

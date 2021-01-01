library(ggplot2)
library(dplyr)    # for data manipulation
library(tidyr) # for pivot
library(data.table)

setwd('~/code/shithouse/prelim_results')

x <- read.csv('studies_consolidated.tsv', sep='\t')
reads.unfiltered <- as.data.frame(x)
# temporarily get rid of the sample names so we can do filtering
rownames(reads.unfiltered) <- reads.unfiltered$X
reads.unfiltered$X <- NULL
reads.unfiltered <- reads.unfiltered[rowSums(reads.unfiltered) > 1000,] # get rid of empty samples
# taxa must have > 70 reads across ALL samples
reads.unfiltered <- reads.unfiltered[,colSums(reads.unfiltered) > 70]
# taxa must appear in a minimum number of samples
#reads = reads[,!sapply(df, function(x) mean(x==0))>0.1]
reads <- reads.unfiltered
reads$srr <- rownames(reads)

colnames(reads) <- gsub('(\\w+\\.\\w+)\\.\\w+\\..+$', '\\1', colnames(reads))
length(unique(colnames(reads)))
reads.longer <- reads %>% pivot_longer(!srr, names_to = "taxon", values_to = "count")

dt <- data.table(reads.longer)
combined <- dt[, lapply(.SD, sum), by=list(srr, taxon)]
final <- combined %>% pivot_wider(names_from=taxon, values_from=count, values_fill=0)
# get rid of sample names, we don't need em anymore
final <- as.data.frame(final)
final$srr <- NULL

rel = function(x, na.rm=TRUE){
  x / sum(x)
}

final.rel <- as.data.frame(t(apply(final, 1, rel)))

rowSums(final.rel) #make sure it's actually adding up

# find the taxa with the highest variance, so we can use them
# to order the samples
variance <- as.data.frame(sapply(final.rel, sd))
variance$taxon <- colnames(final.rel)
colnames(variance) <- c('sd','taxon')
variance <- variance[order(-variance$sd),]

final.rel$sample <- rownames(final.rel)

# if we want to order the samples using certain taxa,
# we need an extra copy of those

# go BACK to long form to plot
final.long <- final.rel %>% pivot_longer(!sample, names_to = "taxon", values_to = "rel")

ggplot(final.long, aes(fill=taxon, y=rel, x=sample)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  theme(
    legend.position='bottom'
  )

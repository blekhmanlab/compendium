library(tidyverse)
library(data.table)

setwd('~/code/shithouse/manuscript')

combine_taxa <- function(dataset) {
  dataset$srr <- rownames(dataset)
  dt <- dataset %>%
    pivot_longer(!srr, names_to = "taxon", values_to = "count") %>%
    data.table()
  final <- dt[, lapply(.SD, sum), by=list(srr, taxon)] %>%
            pivot_wider(names_from=taxon, values_from=count, values_fill=0) %>%
            data.frame() # tibbles don't like row names
  rownames(final) <- final$srr
  final$srr <- NULL
  return(final)
}

get_prevalence <- function(taxtable) {
  richness <- taxtable %>% mutate_if(is.numeric, ~1 * (. > 0))
  prevalence <- data.frame(colnames(richness), colSums(richness))
  colnames(prevalence) <- c('taxon','samples')
  rownames(prevalence) <- NULL
  return(prevalence)
}

rel = function(x){
  x / sum(x)
}

make_rel <- function(table) {
  return(data.frame(t(apply(table, 1, rel))))
}

get_abundance <- function(taxtable) {
  taxtable.rel <- as.data.frame(t(apply(taxtable, 1, rel)))
  abundance <- data.frame(colnames(taxtable.rel), colSums(taxtable.rel, na.rm=TRUE)/nrow(taxtable.rel))
  rownames(abundance) <- NULL
  colnames(abundance) <- c('taxon','mean_abundance')
  return(abundance)
}


#/home/blekhman/shared/compendium/results/taxa_files/

raw <- read.csv('studies_consolidated.tsv', sep='\t')

samplelist <- data.frame(name=gsub('\\w+_consolidated.tsv_(\\w+)$', '\\1', raw$X))
write.csv(samplelist, file='sample_list.csv', row.names = FALSE, quote=FALSE)



x <- raw
rownames(x) <- x$X
x$X <- NULL
nrow(x) # total samples to start
ncol(x) # total taxa to start
nrow(x[rowSums(x) < 1000,]) # how many samples don't have enough reads?
x <- x[rowSums(x) >= 1000,] # get rid of empty samples
rownames(x) <- gsub('(\\w+)_consolidated.tsv(_\\w+)$', '\\1\\2', rownames(x))

# taxa must have > 100 reads across ALL samples
# (this is redundent to the next step, but is useful
# for thinning out the data frame before we start summing
# up columns)
ncol(x[,colSums(x) < 900]) # taxa w less than 100 reads total
x <- x[,colSums(x) >= 900]
# taxa must appear in at least 1000 samples
richness <- x %>% mutate_if(is.numeric, ~1 * (. > 0))
ncol(x[,colSums(richness) < 1000])
x <- x[,colSums(richness) >= 1000]
# after we get rid of taxa, make sure all our sample still have some reads left
nrow(x[rowSums(x) < 1000,])
x <- x[rowSums(x) >= 1000,]
rm(richness)

# Remove the smamples with messed up assignments
taxphylum <- x
colnames(taxphylum) <- gsub('^(\\w+\\.\\w+)\\.\\w+\\..+$', '\\1', colnames(taxphylum))
taxphylum <- combine_taxa(taxphylum) %>% make_rel()
taxphylum$weird <- taxphylum$NA.NA + taxphylum$Eukaryota.NA +
    taxphylum$Bacteria.NA
nrow(taxphylum[taxphylum$weird >= 0.1,]) # samples to remove
x <- x[taxphylum$weird < 0.1,]

nrow(taxphylum[taxphylum$Archaea.Euryarchaeota >= 0.1,]) # samples to remove
x <- x[taxphylum$Archaea.Euryarchaeota < 0.1,]




# project stats
projects <- data.frame(rownames(x)) %>% separate(rownames.x., c("project", "sample"), "_")
projects.count <- data.frame(table(projects$project))
colnames(projects.count) <- c('project','samples')
nrow(projects.count) # how many projects?
median(projects.count$samples)
max(projects.count$samples)

# library size
size <- data.frame(
  sample=rownames(x),
  lib=rowSums(x)
)
mean(size$lib)
sd(size$lib)

check <- x[projects$project=='PRJNA63661',]
colnames(check) <- gsub('^(\\w+\\.\\w+)\\.\\w+\\..+$', '\\1', colnames(check))
check <- check %>%
    combine_taxa() %>%
    make_rel()

# CONSOLIDATE COLUMN NAMES AT DIFFERENT TAXONOMIC LEVELS:
# (todo: remove entries like 'Bacteria.NA.NA.NA.NA' etc)
taxphylum <- x
colnames(taxphylum) <- gsub('^(\\w+\\.\\w+)\\.\\w+\\..+$', '\\1', colnames(taxphylum))
taxphylum <- combine_taxa(taxphylum)
prevalence.phylum <- get_prevalence(taxphylum) # MOST PREVALENT
abundance.phylum <- get_abundance(taxphylum)
# -----------------
taxclass <- x
colnames(taxclass) <- gsub('^(\\w+\\.\\w+\\.\\w+)\\..+$', '\\1', colnames(taxclass))
taxclass <- combine_taxa(taxclass)
prevalence.class<- get_prevalence(taxclass)
abundance.class <- get_abundance(taxclass)
# -----------------
taxorder <- x
colnames(taxorder) <- gsub('^(\\w+\\.\\w+\\.\\w+\\.\\w+)\\..+$', '\\1', colnames(taxorder))
taxorder <- combine_taxa(taxorder)
prevalence.order <- get_prevalence(taxorder)
abundance.order <- get_abundance(taxorder)
# -----------------
taxfamily <- x
colnames(taxfamily) <- gsub('^(\\w+\\.\\w+\\.\\w+\\.\\w+\\.\\w+)\\..+$', '\\1', colnames(taxfamily))
taxfamily <- combine_taxa(taxfamily)
prevalence.family <- get_prevalence(taxfamily)
abundance.family <- get_abundance(taxfamily)
# -----------------
taxgenus <- x
colnames(taxgenus) <- gsub('^(\\w+\\.\\w+\\.\\w+\\.\\w+\\.\\w+.\\w+)', '\\1', colnames(taxgenus))
taxgenus <- combine_taxa(taxgenus)
prevalence.genus <- get_prevalence(taxgenus)
abundance.genus <- get_abundance(taxgenus)

library(tidyverse)
library(data.table)
library(ggrepel)

setwd('~/compendium')

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
x <- read.csv('studies_consolidated.tsv', sep='\t')
rownames(x) <- x$X
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
nrow(taxphylum[taxphylum$weird > 0.1,]) # samples to remove
x <- x[taxphylum$weird <= 0.1,]

# re-calculate this so the rows in taxphylum match the new
# rows in x
taxphylum <- x
colnames(taxphylum) <- gsub('^(\\w+\\.\\w+)\\.\\w+\\..+$', '\\1', colnames(taxphylum))
taxphylum <- combine_taxa(taxphylum) %>% make_rel()
nrow(taxphylum[taxphylum$Archaea.Euryarchaeota > 0.1,]) # samples to remove
x <- x[taxphylum$Archaea.Euryarchaeota <= 0.1,]



# -----------------
taxclass <- x
colnames(taxclass) <- gsub('^(\\w+\\.\\w+\\.\\w+)\\..+$', '\\1', colnames(taxclass))

reads.nozero <- combine_taxa(taxclass)
prevalence.class <- get_prevalence(reads.nozero)
reads.nozero[reads.nozero == 0] <- 1
rm(taxclass)

reads.nozero <- unique(reads.nozero) # GET RID OF DUPLICATES
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
clr <- apply(reads.nozero, 1, function(a) log(a/gm_mean(a)))
clr <- as.data.frame(t(clr))

library(vegan) # for distance calc
library(bigmds) # https://arxiv.org/abs/2007.11919
date()
set.seed(45)
# plotting 8 MDS dimensions:
nmds <- divide_conquer_mds(clr, 10000, 16, 8, n_cores = 4, dist_fn=vegdist, method='euclidean')
date()
saveRDS(nmds, 'nmds_nmds.rds')
points <- data.frame(nmds$points)

colnames(points) <- c('x','y')
# figure out most prevalent taxon in each sample
points$largest <- colnames(clr)[max.col(clr, ties.method = "first")]
points$sample <- rownames(reads.nozero)
colnames(points) <- c('mds1','mds2','mds3','mds4','mds5','mds6','mds7','mds8','largest','sample')
points$project <- gsub('^(\\w+)_\\w+', '\\1', points$sample)

saveRDS(points, 'nmds_points.rds')

# BIPLOT STUFF
correlations <- data.frame(taxon=character(),
                           x=double(),
                           y=double(),
                           stringsAsFactors=FALSE)

# check correlations for the X most prevalence classes:
to_eval <- prevalence.class[order(-prevalence.class$samples),][1:15,]$taxon

for(taxon in to_eval) {
  xcor <- cor(clr[taxon], points$x, method='spearman')
  ycor <- cor(clr[taxon], points$y, method='spearman')
  tresults <- data.frame(
    taxon=taxon,
    x=xcor,
    y=ycor
  )
  correlations <- rbind(correlations, tresults)
}
# just grab the most specific name for the biplot
correlations$taxon <- gsub('^.+\\.(\\w+)$', '\\1', correlations$taxon)

saveRDS(correlations, 'nmds_correlations.rds')

#ggplot(correlations, aes(x=x, y=y,label=taxon,xend=0,yend=0)) +
#  geom_text() +
#  geom_segment() +
#  theme_bw()

correlations.scaled = correlations
correlations.scaled$x <- correlations.scaled$x * 15
correlations.scaled$y <- correlations.scaled$y * 15


pcoa_plot <- ggplot(points, aes(x=x, y=y, color=largest)) +
  geom_point(size=0.1) +
  theme_bw() +
  labs(
    #x=paste('MDS1 (', round(nmds$eigen[1],1), '%)', sep=''), 
    #y=paste('MDS2 (', round(nmds$eigen[2],1), '%)', sep='')
  ) +
  guides(color = guide_legend(override.aes = list(size=10)))


pcoa_plot +
  geom_segment(data=correlations.scaled, aes(x=x, y=y,xend=0,yend=0, color=NULL)) +
  geom_label_repel(data=correlations.scaled,
          aes(x=x, y=y,label=taxon, color=NULL)) +
  theme(legend.position="none")

ggsave('pcoa.pdf')
ggsave('pcoa.png')
#savehistory('./log16jan.txt')


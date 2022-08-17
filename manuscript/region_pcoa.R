library(tidyverse)
library(data.table)
library(ggrepel)
library(tidyr)
library(vegan) # for distance calc
library(bigmds) # https://arxiv.org/abs/2007.11919

setwd('~/code/moved/shithouse/manuscript')

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

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

x <- read.csv('studies_consolidated.tsv', sep='\t')
rownames(x) <- x$X
x$X <- NULL
x <- x[rowSums(x) >= 1000,] # get rid of empty samples
rownames(x) <- gsub('(\\w+)_consolidated.tsv(_\\w+)$', '\\1\\2', rownames(x))

# taxa must have > 900 reads across ALL samples
# (this is redundant to the next step, but is useful
# for thinning out the data frame before we start summing
# up columns)
x <- x[,colSums(x) >= 900]
# taxa must appear in at least 1000 samples
richness <- x %>% mutate_if(is.numeric, ~1 * (. > 0))
x <- x[,colSums(richness) >= 1000]
# after we get rid of taxa, make sure all our sample still have some reads left
x <- x[rowSums(x) >= 1000,]
rm(richness)

# Remove the samples with messed up assignments
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

#reads.nozero$sample <- gsub('.+_(\\w+)$', '\\1', rownames(reads.nozero))
reads.nozero$sample <- rownames(reads.nozero)
countries <- read.csv('sample_countries.csv')
regions <- read.csv('regions.csv')
reads.nozero.countries <- reads.nozero %>%
  left_join(
    countries,
    by='sample'
  ) %>%
  left_join(regions, by=c('country'='alpha2'))

reads.nozero.countries[is.na(reads.nozero.countries$region),]$region <- 'unknown'

# filter to show only samples from region
target <- 'Latin America and the Caribbean'
test <- reads.nozero.countries[reads.nozero.countries$region==target,] %>%
  select(!c('sample','country','region','geo_loc_name','name'))

# we get prevalence now because it's about to get messed up
# by the pseudo-counts
prevalence.class <- get_prevalence(test)
test[test == 0] <- 1

test <- unique(test) # GET RID OF DUPLICATES
clr <- apply(test, 1, function(a) log(a/gm_mean(a)))
clr <- as.data.frame(t(clr))

###### trying PCOA the regular way
#dist <- vegdist(clr, method = "euclidean")
#library(ape)
#PCOA <- pcoa(dist)
#barplot(PCOA$values$Relative_eig[1:10])
#biplot.pcoa(PCOA, clr)
######################

date()
set.seed(45)
# plotting 8 MDS dimensions:
nmds <- divide_conquer_mds(clr, 1000, 16, 8, n_cores = 4, dist_fn=vegdist, method='euclidean')
date()
#saveRDS(nmds, 'nmds_nmds.rds')

points <- data.frame(nmds$points)
# figure out most prevalent taxon in each sample
colnames(points) <- c('mds1','mds2','mds3','mds4','mds5','mds6','mds7','mds8')
#saveRDS(points, 'nmds_points.rds')

# BIPLOT STUFF
correlations <- data.frame(taxon=character(),
                           x=double(),
                           y=double(),
                           stringsAsFactors=FALSE)

# check correlations for the X most prevalent classes:
to_eval <- prevalence.class[order(-prevalence.class$samples),][1:15,]$taxon

for(taxon in to_eval) {
  xcor <- cor(clr[taxon], points$mds1, method='spearman')
  ycor <- cor(clr[taxon], points$mds2, method='spearman')
  tresults <- data.frame(
    taxon=taxon,
    x=xcor,
    y=ycor
  )
  correlations <- rbind(correlations, tresults)
}
# just grab the most specific name for the biplot
correlations$taxon.trimmed <- gsub('^.+\\.(\\w+)$', '\\1', correlations$taxon)

#saveRDS(correlations, 'nmds_correlations.rds')

#ggplot(correlations, aes(x=x, y=y,label=taxon,xend=0,yend=0)) +
#  geom_text() +
#  geom_segment() +
#  theme_bw()

#redoing correlation labels
correlations$taxon <- rownames(correlations)
correlations$taxon <- gsub('^Bacteria.NA.NA$', 'Bacteria (unassigned)', correlations$taxon)
correlations$taxon <- gsub('^Bacteria.Proteobacteria.NA$', 'Proteobacteria', correlations$taxon)
correlations$taxon <- gsub('^NA.NA.NA$', 'Unassigned', correlations$taxon)

correlations$taxon <- gsub('^\\w+\\.\\w+\\.(\\w+)$', '\\1', correlations$taxon)

correlations.scaled <- correlations
correlations.scaled$x <- correlations.scaled$x * 15
correlations.scaled$y <- correlations.scaled$y * 15

scree <- data.frame(
  axis=9-row_number(nmds$eigen),
  eigen=nmds$eigen,
  eigen.rel=100 * (nmds$eigen/sum(nmds$eigen))
)
scree$cumulative <- cumsum(scree$eigen.rel)


ggplot(points, aes(x=mds1, y=mds2) ) +
  geom_hex(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  guides(fill = guide_colorbar(title.position = "top")) +
  labs(
    x=paste('MDS1 (', round(scree$eigen.rel[1],1), '%)', sep=''), 
    y=paste('MDS2 (', round(scree$eigen.rel[2],1), '%)', sep=''), 
    fill="Samples",
    title=target
  ) +
  geom_segment(data=correlations.scaled,
               aes(x=x, y=y,xend=0,yend=0, color=NULL),
               color='white', size=1,
               show.legend = FALSE) +
  geom_label_repel(data=correlations.scaled,
                   aes(x=x, y=y,label=taxon, color=NULL),
                   force_pull=3, show.legend = FALSE)




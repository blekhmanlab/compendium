library(ggplot2)
library(dplyr)    # for data manipulation
library(tidyr) # for pivot
library(data.table)
library(scales) # for y-axis labels
library(ggrepel) # for labels in scatter plot
library(patchwork)

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
rownames(final) <- final$srr
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

# if we want to order the samples using certain taxa,
# we need an extra copy of those, which we'll attach
# to every entry for every sample in the pivot_longer version.
final.rel$sample <- rownames(final.rel)

# now go BACK to long form to plot
final.long <- final.rel %>% pivot_longer(!sample, names_to = "taxon", values_to = "rel")

final.long[!(final.long$taxon %in% variance$taxon[1:5]),]$taxon <- 'other'
final.long$taxon <- factor(final.long$taxon, levels=c(variance$taxon, 'other'))

final.long$sample <- factor(final.long$sample, levels=final.rel[
    order(final.rel[[variance$taxon[1]]],
          final.rel[[variance$taxon[1]]]+final.rel[[variance$taxon[2]]],
          final.rel[[variance$taxon[1]]]+final.rel[[variance$taxon[2]]]+final.rel[[variance$taxon[3]]],
          final.rel[[variance$taxon[1]]]+final.rel[[variance$taxon[2]]]+final.rel[[variance$taxon[3]]]+final.rel[[variance$taxon[4]]],
          final.rel[[variance$taxon[1]]]+final.rel[[variance$taxon[2]]]+final.rel[[variance$taxon[3]]]+final.rel[[variance$taxon[4]]]+final.rel[[variance$taxon[5]]]),
  ]$sample)

# trying to group by more than just the first level makes everything
# look lumpy
#final.long$sample <- factor(final.long$sample, levels=final.rel[
#  order(round(final.rel[[variance$taxon[1]]], digits=1),
#        round(final.rel[[variance$taxon[2]]], digits=1),
#        round(final.rel[[variance$taxon[3]]], digits=1),
#        final.rel[[variance$taxon[4]]]),
#]$sample)

panel_a <- ggplot(final.long, aes(fill=taxon, y=rel, x=sample)) + 
  geom_bar(stat="identity") +
  theme_bw() +
  theme(
    legend.position='none',
    axis.text=element_text(size=11),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  ) +
  scale_fill_brewer(palette='Set1',
     labels=c('Firmicutes', 'Bacteroidota', 'Proteobacteria', 'Actinobacteria','Euryarchaeota','other')) +
  scale_y_continuous(labels = scales::percent, expand=c(0,0)) +
  labs(x='Sample',y='Relative abundance', fill='Class')

a_legend <- cowplot::get_legend(panel_a +
              theme(
                legend.position='bottom',
                legend.text = element_text(size=11)
              )
            )
plot(a_legend)


totalcounts <- data.frame(
  continent = c('North America','Europe', 'Asia','Africa','Oceania','South America','(unknown)'),
  samples = c(78123,42939,22201,9479,4636,3843,8797),
  label_pos = c(71123,50939,30201,15479,10636,9843,14797)
)
totalcounts$continent <- factor(totalcounts$continent,
                                   levels=rev(totalcounts$continent))

panel_b <- ggplot(totalcounts, aes(x=continent, y=samples, fill=continent)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=c("#AEAEAE","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#00CFFF")) +
  geom_text(aes(x=continent, y=label_pos,
      label=scales::comma(samples),
      color=c(rep('black',7))),
      fontface='bold') +
    scale_color_manual(values=c('black','white')) +
  theme_bw() +
  labs(x='Continent', y='Human fecal samples') +
  coord_flip() +
  scale_y_continuous(labels=comma, limits=c(0,85000), expand=c(0, 0)) +
  theme(
    axis.text=element_text(size=11),
    legend.position='none'
  )

# t-SNE
library(Rtsne)
load(file="reads.nozero.Rda")
reads.nozero <- unique(reads.nozero) # GET RID OF DUPLICATES

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
clr <- apply(reads.nozero, 1, function(a) log(a/gm_mean(a)))
clr <- as.data.frame(t(clr))

#to_tsne$sample <- NULL
#to_tsne <- unique(to_tsne)

set.seed(1234)
tsned <- Rtsne(as.matrix(clr), perplexity=475,
               max_iter=5000, num_threads=3,
               theta=0.5)
toplot <- as.data.frame(tsned$Y)

allcountries <- read.csv('sample_countries.csv')

countries_toplot <- as.data.frame(rownames(clr))
colnames(countries_toplot) <- c('srr')
# trim off project names
countries_toplot$srr <- gsub('^[^_]+_','', countries_toplot$srr)
countries_toplot <- countries_toplot %>% left_join(allcountries, by = "srr")
countries_toplot[is.na(countries_toplot$X),]$X <- '?'
countries_toplot[countries_toplot$X=='?',]$X <- '(unknown)'

panel_c <- ggplot(toplot, aes(x=V1, y=V2, color=countries_toplot$X)) +
  geom_point() +
  #scale_color_manual(values=c("#AEAEAE","#D95F02","#7570B3","#66A61E","#E6AB02","#00CFFF")) +
  scale_color_manual(values=c("#AEAEAE","#66A61E","#E6AB02","#00CFFF","#7570B3","#D95F02")) +
  labs(x='tSNE1',y='tSNE2',color='Continent') +
  theme_bw() +
  theme(
    legend.position='none',
    axis.text=element_text(size=11)
  )


panel_a / (panel_b|panel_c) +
  plot_annotation(tag_levels = 'A')


#TODO: ADD CONTINENT FACET??

# prevalence
prevalence <- as.data.frame(colSums(final.rel != 0)/nrow(final.rel))
colnames(prevalence) <- c('prevalence')

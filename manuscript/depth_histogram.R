library(ggplot2)

setwd('~/code/shithouse/manuscript')

x <- read.csv('studies_consolidated.tsv', sep='\t')
unfiltered <- as.data.frame(x)
unfiltered$X <- NULL

unfiltered$depth <- rowSums(unfiltered)

ggplot(unfiltered, aes(x=depth)) + 
  geom_histogram(binwidth=10000) +
  theme_bw() +
  scale_y_log10()

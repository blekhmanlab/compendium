library(ggplot2)
library(scales) # for y-axis labels
library(ggrepel) # for labels in scatter plot
library(patchwork)

final.rel <- make_rel(taxphylum)

# find the taxa with the highest variance, so we can use them
# to order the samples
variance <- as.data.frame(sapply(final.rel, sd))
variance$taxon <- colnames(final.rel)
colnames(variance) <- c('sd','taxon')
variance <- variance[order(-variance$sd),]

# we can't possibly display 180,000 samples on a screen
# that's MAYBE 3000 pixels wide, so grab a random sample
# of 5,000:
set.seed(42)
subset <- final.rel[sample(nrow(final.rel), 5000), ]

# if we want to order the samples using certain taxa,
# we need an extra copy of those, which we'll attach
# to every entry for every sample in the pivot_longer version.
subset$sample <- rownames(subset)

# now go BACK to long form to plot
final.long <- subset %>% pivot_longer(!sample, names_to = "taxon", values_to = "rel")

final.long[!(final.long$taxon %in% variance$taxon[1:5]),]$taxon <- 'other'
final.long$taxon <- factor(final.long$taxon, levels=c(variance$taxon, 'other'))

final.long$sample <- factor(final.long$sample, levels=subset[
    order(subset[[variance$taxon[1]]],
          subset[[variance$taxon[1]]]+subset[[variance$taxon[2]]],
          subset[[variance$taxon[1]]]+subset[[variance$taxon[2]]]+subset[[variance$taxon[3]]],
          subset[[variance$taxon[1]]]+subset[[variance$taxon[2]]]+subset[[variance$taxon[3]]]+subset[[variance$taxon[4]]],
          subset[[variance$taxon[1]]]+subset[[variance$taxon[2]]]+subset[[variance$taxon[3]]]+subset[[variance$taxon[4]]]+subset[[variance$taxon[5]]]),
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
  geom_bar(stat="identity",width=1) +
  theme_bw() +
  theme(
    legend.position='none',
    axis.text=element_text(size=11),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  ) +
  scale_fill_brewer(palette='Set1',
     labels=c('Firmicutes', 'Bacteroidota', 'Proteobacteria', 'Actinobacteria','Verrucomicrobiota','other')) +
  scale_y_continuous(labels = scales::percent, expand=c(0,0)) +
  labs(x='Sample',y='Relative abundance', fill='Phylum')

a_legend <- cowplot::get_legend(panel_a +
              theme(
                legend.position='bottom',
                legend.text = element_text(size=11),
                plot.tag.position  = c(1, 1)
              )
            )

#------------- SHANNON DIVERSITY
library(vegan)
library(dplyr)

shannon <- diversity(taxfamily, "shannon")
shannon <- as.data.frame(shannon)
shannon$srr <- rownames(shannon)

depth <- data.frame(rownames(taxfamily), rowSums(taxfamily))
colnames(depth) <- c('srr','depth')

diversity <- ggplot(shannon, aes(x=shannon)) + 
  geom_histogram(bins=50, fill='#969696', color='black') +
  geom_vline(xintercept=median(shannon$shannon),
             size=1,color="red") +
  annotate("text", label=paste(
      "median: ", round(median(shannon$shannon), 3)
    ),
    size = 4,
    x=1.4, y=9000
  ) +
  theme_bw() +
  labs(
    x="Shannon diversity",
    y="Samples"
  ) +
  scale_y_continuous(labels=scales::comma)

depth_plot <- ggplot(depth, aes(x=depth)) + 
  geom_histogram(bins=50, fill='#969696', color='black') +
  geom_vline(xintercept=median(depth$depth),
             size=1,color="red") +
  annotate("text", label=paste(
      "median: ", median(depth$depth)
    ),
    size = 4,
    x=80000, y=14000
  ) +
  theme_bw() +
  labs(
    x="Merged reads (log)",
    y="Samples"
  ) +
  scale_y_continuous(
    labels=scales::comma,
    limits=c(0,15000),
    expand=c(0,0)
  ) +
  scale_x_log10(labels=scales::comma)

#------- assembly
big_stacked <- panel_a 

(panel_a + inset_element(a_legend, align_to='plot',
            0.6, 0.6, 0.95, 0.8,
            ignore_tag = TRUE) ) / (diversity | depth_plot) +
  plot_annotation(tag_levels='A')


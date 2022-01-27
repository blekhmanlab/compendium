library(ggplot2)
library(scales) # for y-axis labels
library(ggrepel) # for labels in scatter plot
library(patchwork)

final.rel <- as.data.frame(t(apply(taxphylum, 1, rel)))

rowSums(final.rel) #make sure it's actually adding up

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
                legend.text = element_text(size=11)
              )
            )

panel_a + inset_element(a_legend, 0.6, 0.8, 0.9, 0.8)

library(ggplot2)
library(ggrepel)
library(ggExtra) #for histograms

nmds <- readRDS('nmds_nmds.rds')
points <- readRDS('nmds_points.rds')
correlations <- readRDS('nmds_correlations.rds')

#redoing correlation labels
correlations$taxon <- rownames(correlations)
correlations$taxon <- gsub('^Bacteria.NA.NA$', 'Unassigned', correlations$taxon)
correlations$taxon <- gsub('^Bacteria.Firmicutes.NA$', 'Firmicutes', correlations$taxon)
correlations$taxon <- gsub('^\\w+\\.\\w+\\.(\\w+)$', '\\1', correlations$taxon)

correlations.scaled <- correlations
correlations.scaled$x <- correlations.scaled$x * 12
correlations.scaled$y <- correlations.scaled$y * 12

correlations.scaled$dist <- sqrt(abs(correlations.scaled$x) + abs(correlations.scaled$y))

correlations.scaled <- correlations.scaled %>%
  arrange(desc(dist)) %>%
  slice_head(n=10)

# THIS IS FOR DEBUGGING WEIRD SECTIONS OF THE DIAGRAM:
#trimpoints <- points

#trimpoints <- points[points$mds1 > 10,]
#trimpoints <- trimpoints[trimpoints$mds2 < -5.5,]

#ggplot(trimpoints, aes(x=mds1, y=mds2, color=largest)) +
#  geom_point(size=0.1) +
#  theme_bw() +
#  guides(color = guide_legend(override.aes = list(size=10)))

#answer <- data.frame(table(trimpoints$project))

# Make the labels pretty, chop off "Bacteria" at the beginning of each
points$largest <- gsub('^(\\w+)\\.(\\w+)\\.(\\w+)$', '\\2', points$largest)

points$largest <- gsub('^\\w+\\.(\\w+)\\.(\\w+)$', '\\1 \\2', points$largest)
points$largest <- gsub(' NA$', ' (Unassigned)', points$largest)

left <- ggplot(points, aes(x=x, y=y, color=largest)) +
  geom_point(size=0.1) +
  theme_bw() +
  labs(
    x=paste('MDS1 (', round(nmds$eigen[1],1), '%)', sep=''), 
    y=paste('MDS2 (', round(nmds$eigen[2],1), '%)', sep=''),
    color="Most abundant phylum"
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size=10),
      title.position = "top"
    )
  ) +
  theme(legend.position="bottom") +
  geom_segment(data=correlations.scaled, aes(x=x, y=y,xend=0,yend=0, color=NULL),
               show.legend = FALSE) +
  geom_label_repel(data=correlations.scaled,
                   aes(x=x, y=y,label=taxon, color=NULL),
                   force_pull=3, show.legend = FALSE)


right <- ggplot(points, aes(x=x, y=y) ) +
  geom_hex(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() +
  guides(fill = guide_colorbar(title.position = "top")) +
  labs(
    x=paste('MDS1 (', round(nmds$eigen[1],1), '%)', sep=''), 
    y=element_blank(),
    fill="Samples"
  ) +
  geom_segment(data=correlations.scaled,
               aes(x=x, y=y,xend=0,yend=0, color=NULL),
               color='white', size=1,
               show.legend = FALSE) +
  geom_label_repel(data=correlations.scaled,
                   aes(x=x, y=y,label=taxon, color=NULL),
                  force_pull=3, show.legend = FALSE) +
  theme(legend.position='bottom')

layout <- "
AAABBB
CC##DD
"

(left | right)  +
  plot_layout(guides='collect') +
  plot_annotation(tag_levels = 'A') &
  theme(legend.position='bottom')

ggsave('pcoa.pdf',
       plot=p2,
       width=3000,
       height=5000,
       units="px",
       device="pdf"
)

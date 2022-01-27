library(ggplot2)
library(ggrepel)

nmds <- readRDS('nmds_nmds.rds')
points <- readRDS('nmds_points.rds')
correlations <- readRDS('nmds_correlations.rds')

correlations.scaled <- correlations
correlations.scaled$x <- correlations.scaled$x * 15
correlations.scaled$y <- correlations.scaled$y * 15


# THIS IS FOR DEBUGGING WEIRD SECTIONS OF THE DIAGRAM:
#trimpoints <- points

#trimpoints <- points[points$mds1 > 10,]
#trimpoints <- trimpoints[trimpoints$mds2 < -5.5,]

#ggplot(trimpoints, aes(x=mds1, y=mds2, color=largest)) +
#  geom_point(size=0.1) +
#  theme_bw() +
#  guides(color = guide_legend(override.aes = list(size=10)))

#answer <- data.frame(table(trimpoints$project))

pcoa_plot <- ggplot(points, aes(x=mds1, y=mds2, color=largest)) +
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
                   aes(x=x, y=y,label=taxon, color=NULL),
                  force_pull=3) +
  theme(legend.position="none")

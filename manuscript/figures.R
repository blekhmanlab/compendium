library(cowplot)

(taxa_bars | pcoa_fig ) / big_stacked


FIG_ONE <- plot_grid(
  plot_grid(taxa_bars, pcoa_fig,nrow = 1),
  big_stacked,
  ncol=1, rel_heights=c(3,2)
)

FIG_ONE
ggsave('fig1.png',
    plot=FIG_ONE,
    width=2250,
    height=2500,
    units="px",
    device="png"
)

# ------------------------------------------------------
##
## ---- ggplot theme ----
##
# ------------------------------------------------------
## This script should run fine if you have run the
## Functions_Source-for-CDBN-Analyses R script first.
require(tidyverse)
require(gridExtra)

## Many thanks to Rob Heckman for his ggplot theme script!

theme_oeco <- theme_classic() +
  theme(axis.title = element_text(size = 10), axis.text = element_text(size = 10),
        axis.line.x = element_line(size = 0.35, colour = 'grey50'),
        axis.line.y = element_line(size = 0.35, colour = 'grey50'),
        axis.ticks = element_line(size = 0.25, colour = 'grey50'),
        legend.justification = c(1, 0.75), legend.position = c(1, 0.9),
        legend.key.size = unit(0.35, 'cm'),
        legend.title = element_blank(), legend.text = element_text(size = 9),
        legend.text.align = 0, legend.background = element_blank(),
        plot.subtitle = element_text(size = 10, vjust = 0), #plot.margin = unit(c(0.35, 0, 0.25, 0), 'cm'),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0.5, size = 10 ,vjust = 0),
        strip.placement = 'outside', panel.spacing.x = unit(-0.5, 'cm'))
  theme_set(theme_oeco)

theme_poster <- theme_classic() +
    theme(axis.title = element_text(size = 18), axis.text = element_text(size = 14),
          axis.line.x = element_line(size = 0.35, colour = 'grey50'), axis.line.y = element_line(size = 0.35, colour = 'grey50'),
          axis.ticks = element_line(size = 0.5, colour = 'grey50'),
          legend.justification = c(1, 0.75), legend.position = c(1, 0.9), legend.key.size = unit(0.5, 'cm'),
          legend.title = element_blank(), legend.text = element_text(size = 18),
          legend.text.align = 0, legend.background = element_blank(),
          plot.subtitle = element_text(size = 14, vjust = 0), #plot.margin = unit(c(0.35, 0, 0.25, 0), 'cm'),
          strip.background = element_blank(), strip.text = element_text(hjust = 0.5, size = 10 ,vjust = 0),
          strip.placement = 'outside', panel.spacing.x = unit(-0.5, 'cm'))
#  theme_set(theme_poster)


update_geom_defaults('boxplot', list(size = 0.5))
update_geom_defaults('pointrange', list(size = 0.35))
update_geom_defaults('smooth', list(span = 1, size = 0.5))
update_geom_defaults('point', list(size = 2))
update_geom_defaults('errorbar', list(size = 0.35))
update_geom_defaults('hline', list(size = 0.35, colour = 'grey50'))

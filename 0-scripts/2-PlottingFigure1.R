
### load packages
rm(list = ls())
library(tidyverse)
library(gghalves)
library(cowplot)

phenotype = read.csv('../1-PhenotypicAnalysis&GenstatFormatting/Ionome_cleaned.csv')
h2 = read.csv('..//2-Heritability&GeneticCorrelation/Heritability.csv')

phenotype2 = phenotype %>% gather('Traits','Values',-SITE,-LINE) %>% 
  mutate(Type = case_when(
    str_detect(Traits,'Ca|Mg|P|K') ~ 'macronutrient',
    str_detect(Traits,'Sr|Rb') ~ 'analogue',
    str_detect(Traits,'Mn|Zn|Cu|Mo|Fe|B|Co|Se')~'micronutrient',
    str_detect(Traits,'Al|Na|As|Cd') ~ 'harmful'
  ))

theme_Li = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.title = element_text(hjust = 0.5),
  axis.title = element_text(size = 10), axis.text = element_text(size = 8),
  axis.line.x = element_line(size = 0.35, colour = 'grey50'),
  axis.line.y = element_line(size = 0.35, colour = 'grey50'),
  axis.ticks = element_line(size = 0.25, colour = 'grey50'),
  strip.background = element_blank(),
  strip.text = element_text(hjust = 0.5, size = 8 ,vjust = 0)
)

phenotype2$Type = factor(phenotype2$Type, levels = c('macronutrient','analogue','micronutrient','harmful'))

p1 =  ggplot(phenotype2, aes(x= SITE, y=Values, fill=SITE))+geom_half_violin()+
  facet_wrap(~Type+Traits, nrow=6, scales = 'free')+theme_bw()+
  theme_Li+ xlab("")+ylab("") + theme(legend.position="none")
  
h2 = h2 %>% mutate(Type = case_when(
  str_detect(TRAIT,'Ca|Mg|P|K') ~ 'macronutrient',
  str_detect(TRAIT,'Sr|Rb') ~ 'analogue',
  str_detect(TRAIT,'Mn|Zn|Cu|Mo|Fe|B|Co|Se')~'micronutrient',
  str_detect(TRAIT,'Al|Na|As|Cd') ~ 'harmful'
))
h2$Type = factor(h2$Type, levels = c('macronutrient','analogue','micronutrient','harmful'))

p2 = ggplot(h2,aes(SITE,h2))+geom_col(alpha=0.6)+ facet_wrap(~Type+TRAIT, nrow=6,scales = 'free')+
  theme_bw()+ 
  theme_Li+
  ylab('Heritability')+xlab('')+
  geom_errorbar(aes(ymin=h2-h2_SE, ymax=h2+h2_SE), width=0)

F1 = plot_grid(p1,p2, labels = c('(a)','(b)'),nrow = 2)
save_plot(file='../1-PhenotypicAnalysis&GenstatFormatting/Fig1_PhenotypicVariation&heritability.pdf',
          F1, base_height = 18, base_width = 8)



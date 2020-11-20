### Formatting Genstat Results##
rm(list=ls())
library(tidyverse)
library(reshape2)
library(xlsx)
library(ape)
library(ComplexUpset)
library(cowplot)
setwd("../2-GenstatRunning/")

###1.formatting EFF files, get the QTL effects
QTLEffectFiles=intersect(list.files(pattern='EFF'),list.files(pattern='xlsx'))
qtlEff = lapply(QTLEffectFiles, function(x){
  tmp = read.xlsx(x, sheetIndex = 1)
  tmp2 = tmp %>% gather('ENVI','VALUE',-POS) %>% separate(ENVI, c('TYPE','CROSS','SITE')) %>% dplyr::rename(INDEX=POS)%>%
    spread(TYPE,VALUE) %>% mutate(TRAIT = strsplit(strsplit(x,"_",fixed=T)[[1]][2],".", fixed= T)[[1]][1] )
})
qtlEff = do.call(rbind, qtlEff)

qtlEff = qtlEff %>% mutate(SITE=str_replace(SITE,'PKLE','TX'))%>% 
  mutate(SITE=str_replace(SITE,'CLMB','MO'))%>% 
  mutate(SITE=str_replace(SITE,'KBSM','MI'))%>%
  mutate(TRAIT=str_replace(TRAIT,'[0-9]+',""))

write.csv(qtlEff,'../3-GenstatOutput/QTLEffect.csv', row.names = F)

###2. pairwised t test for the significance of 
####differential sensitivity (DS) and atagonistic pleiotropy (AP)
eff2 = qtlEff %>% filter(complete.cases(.))

Eff_Sig = NULL
for (phe in unique(eff2$TRAIT)){
  tmp1 = eff2 %>% filter(TRAIT==phe)
  for (i in unique(tmp1$INDEX)){
    tmp2 = tmp1 %>% filter(INDEX==i)
    for (cross in unique(tmp2$CROSS)){
      tmp3 = tmp2 %>% filter(CROSS==cross)
      
      Z_mo_mi = (tmp3$EFF[1]-tmp3$EFF[2])/(sqrt(tmp3$SE[1]^2 + tmp3$SE[2]^2))
      P_mo_mi = (1-pnorm(abs(Z_mo_mi)))*2
      
      Z_mo_tx = (tmp3$EFF[1]-tmp3$EFF[3])/(sqrt(tmp3$SE[1]^2 + tmp3$SE[3]^2))
      P_mo_tx = (1-pnorm(abs(Z_mo_tx)))*2
      
      Z_mi_tx = (tmp3$EFF[2]-tmp3$EFF[3])/(sqrt(tmp3$SE[2]^2 + tmp3$SE[3]^2))
      P_mi_tx = (1-pnorm(abs(Z_mi_tx)))*2
      
      Pvalue = data.frame(Pvalue =c(P_mo_mi,P_mo_tx,P_mi_tx))
      tmp3 = cbind(tmp3, Pvalue)
      
      if (tmp3$Pvalue[1] < 0.05 |tmp3$Pvalue[2]<0.05 |tmp3$Pvalue[3] <0.05 )
      {Sig = data.frame(Sig = c('Yes','Yes','Yes'))}else{Sig = data.frame(Sig = c('No','No','No'))}
      tmp3 = cbind(tmp3, Sig)
      
      if (tmp3$EFF[1] == tmp3$EFF[2] ){GxE = data.frame(GxE = c('NO','NO','NO'))
      } else if(tmp3$EFF[1] >0 & tmp3$EFF[2] >0 & tmp3$EFF[3] >0 & tmp3$EFF[1]!= tmp3$EFF[2]){
        GxE = data.frame(GxE = c('DS','DS','DS'))
      } else if(tmp3$EFF[1] <0 & tmp3$EFF[2] <0 & tmp3$EFF[3] <0 & tmp3$EFF[1]!= tmp3$EFF[2] ){
        GxE = data.frame(GxE = c('DS','DS','DS'))
      } else {GxE = data.frame(GxE = c('AP','AP','AP'))}
      
      tmp3 = cbind(tmp3, GxE)
      Eff_Sig = rbind(Eff_Sig,tmp3)
    }
  }
}

write.csv(Eff_Sig,'Significance_AP_DS.csv', row.names = F)

###3.formating LOD files
genomap = read.delim('4way_male_reduced_marker.map',sep = '', header = F)
genomap = genomap %>% filter(!str_detect(V1,'group')) %>% dplyr::rename(MARKER=V1, POS = V2) %>% 
  mutate_at(vars(POS, MARKER),list(as.character)) %>% mutate_at(vars(POS),list(as.numeric))

qtlLod = lapply(intersect(list.files(pattern = "LOD"), list.files(pattern='.xlsx')), function(x){
  tmp = read.xlsx(x, sheetIndex = 1)
  tmp2 = tmp %>% mutate(MARKER = genomap$MARKER, POS=genomap$POS, TRAIT = strsplit(strsplit(x,"_",fixed=T)[[1]][2],".", fixed= T)[[1]][1] ) %>%
    dplyr::rename(LOD= qtl_minlogp_) %>% dplyr::select(MARKER, POS, LOD, TRAIT) %>% rownames_to_column('INDEX') %>% mutate_at(vars(INDEX),list(as.numeric))
})
qtlLod = do.call(rbind, qtlLod)
qtlLod = qtlLod %>% mutate(TRAIT=str_replace(TRAIT,'[0-9]+',""))

QTLlist = qtlEff %>% left_join(qtlLod) %>% dplyr::select(INDEX,TRAIT, MARKER,POS, LOD) %>% unique() %>% 
  mutate_at(vars(POS, MARKER),list(as.character))  %>% rowwise() %>% mutate(CHR = substr(MARKER, 5,6))%>%
  mutate_at(vars(POS),list(as.numeric)) %>% mutate(QTL = paste0(CHR,'@', round(POS,2))) 

###4. add the flank marker for these QTLs
QTL_flank = NULL
for (phe in unique(QTLlist$TRAIT)){
  tmp1 = QTLlist %>% filter(TRAIT==phe) 
  tmp2 = qtlLod %>% filter(TRAIT==phe) 
  
  for (i in 1:nrow(tmp1)){
    tmp3 = tmp2 %>% filter(MARKER== as.character(tmp1[i,'MARKER']))
    threshold = tmp3$LOD - 1.5
    m = tmp1$INDEX[i] -1 ##one marker lower than the marker with the highest LOD
    n = tmp1$INDEX[i] +1 ##one marker higher than the marker with the highest LOD
    
    while(tmp2$LOD[m] > threshold) {m=m-1}
    tmp1$lowposition[i] = tmp2$POS[m]
    tmp1$down_marker[i] = tmp2$MARKER[m]
    
    while(tmp2$LOD[n] > threshold){n=n+1}
    tmp1$highposition[i] = tmp2$POS[n]
    tmp1$up_marker[i] = tmp2$MARKER[n]
    
    if (tmp1$lowposition[i]> tmp1$POS[i]){
      tmp1$lowposition[i] = tmp2$POS[m+1]
      tmp1$down_marker[i] = tmp2$MARKER[m+1]
    }
    
    if (tmp1$highposition[i]==0){
      tmp1$highposition[i] = tmp2$POS[n-1]
      tmp1$up_marker[i] = tmp2$MARKER[n-1]
    }
    
    tmp4 = tmp2 %>% filter(str_detect(MARKER, tmp1$CHR[i]))
    if (tmp1$highposition[i]< tmp1$lowposition[i] & tmp1$POS[i] > 30){
      tmp1$highposition[i] = tmp4$POS[nrow(tmp4)]
      tmp1$up_marker[i] = tmp4$MARKER[nrow(tmp4)]
    }
   
    if (tmp1$highposition[i]< tmp1$lowposition[i] & tmp1$POS[i] < 30){
      tmp1$lowposition[i] = tmp4$POS[1]
      tmp1$down_marker[i] = tmp4$MARKER[1]
    }
  }
  
  QTL_flank = rbind(QTL_flank,tmp1)
}

###5. check if there is g x e, using qtlEFF file
for (i in 1:nrow(QTLlist)){
  tmp = qtlEff %>% filter(INDEX==QTLlist$INDEX[i] & TRAIT==QTLlist$TRAIT[i] & CROSS=='add')
  if (tmp$EFF[1] != tmp$EFF[2]){QTL_flank$QxE[i] = 'Yes'}
  if (tmp$EFF[1] == tmp$EFF[2]){QTL_flank$QxE[i] = 'No'}
}

QTL_flank = QTL_flank %>% mutate(ChrNo = str_extract(CHR, '[1-9]'), ChrAnnot = str_extract(CHR, '[A-Z]')) %>%
  mutate_at(vars('ChrNo'), list(as.numeric)) %>%  mutate(Dist = case_when(
    str_detect(ChrAnnot,'K') ~ (ChrNo*2-1),
    str_detect(ChrAnnot,'N') ~ ChrNo*2
  ))
write.csv(QTL_flank,'../3-GenstatOutput/QTLlistForAllTraits_flankmarker.csv',row.names = F)


###6. plot summary of QTLs and Upset Plot for QTL colocalization
load('../0-data/sexAveraged.cross.rda')
QTL_flank$TRAIT = factor(QTL_flank$TRAIT)
cols  =  rainbow(n=14)
library(qtlTools)
pdf('Fig2_1_QTL_summary.pdf', width = 14)
with(QTL_flank, segmentsOnMap(cross, phe = TRAIT, chr = CHR, 
                              l = lowposition, h = highposition,
                              peaklod = LOD, peakcM = POS,  
                              showPeaks = TRUE,
                              chrBuffer = c(0.15,0.1), 
                              tick.width=0.05, lwd="byLod",
                              col = cols,
                              leg.inset=-0.006, legendCex=0.7, 
                              legendPosition="right"))

for (phe in unique(QTL_flank$TRAIT)){
  tmp = subset(QTL_flank, TRAIT==phe)
  for (i in 1:nrow(tmp)){
    if (tmp$QxE[i]=='Yes')
      text(x=(tmp$Dist[i]+0.4), y=(tmp$lowposition[i]-2),'*',col = 'red',cex = 1.5)
  }
}
dev.off()

###7. did some manual work here to figure out the overlaps between elements
qtldf2 = read.csv('QTLlist_for_overlaps.csv')
qtl_v <- c("Mg", "P", "K", "Ca", "Rb", "Sr", "Mn", "Fe", "Cu", "Zn", "Mo",
           "Na", "Al", "Cd")
qtl_u =upset(data = qtldf2, intersect = qtl_v, 
             name="Element QTL Colocalization by QTL",
             themes=upset_modify_themes(
               list(
                 'overall_sizes'=theme(axis.text.x=element_text(angle=45, hjust =1))
               ))) 
save_plot(filename = "Fig2_2_Upset_plot_Element_QTL_colocalization.pdf", plot = qtl_u, base_height = 6)

###put these two figures together (manually)---Figure 2 for manuscript

###8. plotting QTL effects and SE
qtlEff2 = qtlEff %>% filter(CROSS != 'dom')%>% mutate(CROSS = str_replace(CROSS,'add2','C x D')) %>% 
  mutate(CROSS = str_replace(CROSS,'add','A x B'))  %>% left_join(QTLlist) 

qtlEff2$SITE = factor(qtlEff2$SITE, levels = c('TX','MO','MI'))

###plot for Na, Mn, Rb, and P only for manuscript Figure 3
## define the theme for plots
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

myplots_Fig3_1 = vector('list',0)
for (phe in c('Na','Mn','Rb')){
  tmp1 = subset(qtlEff2, TRAIT==phe)
  myplots_Fig2_1[[phe]] =
    ggplot(tmp1, aes(SITE, EFF)) + geom_col() + facet_grid(CROSS~QTL) + ggtitle(phe)+
    geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
    theme_Li+
    geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), size=0.5, width=0) + 
    ylab('QTL Effect') +xlab("")
}
F3_1=plot_grid(plotlist =  myplots_Fig3_1, nrow=1, rel_widths = c(0.4,0.8,1.2))

myplots_Fig3_2 = vector('list',0)
for (phe in c('P')){
  tmp1 = subset(qtlEff2, TRAIT==phe)
  myplots_Fig2_2[[phe]] =
    ggplot(tmp1, aes(SITE, EFF)) + geom_col() + facet_grid(CROSS~QTL) + ggtitle(phe)+
    geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
    theme_Li+
    geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), size=0.5, width=0) + 
    ylab('QTL Effect') +xlab("")
}
F3_2 = plot_grid(plotlist =  myplots_Fig3_2, nrow=1)
F3 = gridExtra::grid.arrange(F3_1,F3_2, nrow=2)

ggsave(file='../3-GenstatOutput/Fig3_QTLeffectPlot.pdf',p3, height = 6, width = 12)

###plot all elements for Supplemental Figure S1, macronutrients, 
###analogues, micronutrients, and harmful

myplot1 = vector('list',0)
for (phe in c('Ca','P')){
  tmp1 = subset(qtlEff2, TRAIT==phe)
  myplot1[[phe]] =
    ggplot(tmp1, aes(SITE, EFF)) + geom_col() + facet_grid(CROSS~QTL) + ggtitle(phe)+
    geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
    theme_Li+
    geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), size=0.5, width=0) + 
    ylab('QTL Effect') +xlab("")
}
p1 = plot_grid(plotlist =  myplot1, nrow=1, rel_widths = c(0.3,1))

myplot2 = vector('list',0)
for (phe in c('K','Mg')){
  tmp1 = subset(qtlEff2, TRAIT==phe)
  myplot2[[phe]] =
    ggplot(tmp1, aes(SITE, EFF)) + geom_col() + facet_grid(CROSS~QTL) + ggtitle(phe)+
    geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
    theme_Li+
    geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), size=0.5, width=0) + 
    ylab('QTL Effect') +xlab("")
}
p2 = plot_grid(plotlist =  myplot2, nrow=1, rel_widths = c(0.5,0.8))

myplot3 = vector('list',0)
for (phe in c('Sr','Rb')){
  tmp1 = subset(qtlEff2, TRAIT==phe)
  myplot3[[phe]] =
    ggplot(tmp1, aes(SITE, EFF)) + geom_col() + facet_grid(CROSS~QTL) + ggtitle(phe)+
    geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
    theme_Li+
    geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), size=0.5, width=0) + 
    ylab('QTL Effect') +xlab("")
}
p3=plot_grid(plotlist =  myplot3, nrow=1, rel_widths = c(0.5,0.55))

myplot4 = vector('list',0)
for (phe in c('Fe','Mo','Zn','Mn','Cu')){
  tmp1 = subset(qtlEff2, TRAIT==phe)
  myplot4[[phe]] =
    ggplot(tmp1, aes(SITE, EFF)) + geom_col() + facet_grid(CROSS~QTL) + ggtitle(phe)+
    geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
    theme_Li+
    geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), size=0.5, width=0) + 
    ylab('QTL Effect') +xlab("")
}
p4 = plot_grid(plotlist =  myplot4, nrow=1, rel_widths = c(0.15,0.15,0.2,0.3,0.3))

myplot5 = vector('list',0)
for (phe in c('Na','Cd','Al')){
  tmp1 = subset(qtlEff2, TRAIT==phe)
  myplot5[[phe]] =
    ggplot(tmp1, aes(SITE, EFF)) + geom_col() + facet_grid(CROSS~QTL) + ggtitle(phe)+
    geom_hline(yintercept = 0,linetype=2)+ xlab("")+ theme_bw()+
    theme_Li+
    geom_errorbar(aes(ymin=EFF-SE, ymax=EFF+SE), size=0.5, width=0) + 
    ylab('QTL Effect') +xlab("")
}
p5 = plot_grid(plotlist =  myplot5, nrow=1, rel_widths = c(1,1,2))

FS1 = gridExtra::grid.arrange(p1,p2,p3,p4,p5, nrow=5)

ggsave(file='../3-GenstatOutput/Supplemental_Fig_S1_QTLeffectPlot.pdf',FS1, height = 8, width = 14)


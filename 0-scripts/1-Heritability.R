
### load packages
rm(list = ls())
library(tidyverse)
library(data.table)
library(pegas)
library(xlsx)
library(sommer)

###0. set working directory, should have genotyping file and phenotyping file
setwd(dir = 'c://Users/Li Zhang/Desktop/fourway-ionomics/0-data/')
getwd()

###1. read the genotype .qua file to get the genotype ID
genotype = read.loci('4way.qua',skip = 4)
genotype = genotype %>% dplyr::select(id) %>% dplyr::rename(LINE=id)

###2. read the raw phenotype file
phenotype_raw = read.csv('Ionome.csv')
plantlist = readxl::read_excel('NSF_4WCR_Master Plant List_Final.xlsx')

###removing negative values is traits-based,  set as NA
##3. formating to genstat and get the mean value for each line
phenotype_avg = phenotype_raw %>% left_join(plantlist, by=c('sample'='PLOT_GL')) %>%
  dplyr::select(SITE, LINE,SampleWeight:Cd111) %>% mutate_all(., list(~ifelse(.<0, NA, .)))%>% 
  group_by(SITE, LINE) %>% summarise_all(list(mean),na.rm=T) %>% ungroup()

phenotype_ordered = vector('list', 0)
for (s in unique(phenotype_avg$SITE)){
  phenotype_ordered[[s]] = phenotype_avg %>%  filter(SITE==s) %>% right_join(genotype)%>% mutate(SITE = s)
}
phenotype_ordered =bind_rows(phenotype_ordered)

###load files needed for running sommer
load('c:/Users/Li Zhang/Desktop/fourway-ionomics/0-data/DataFilesNeededForSommer.Rdata')

########heritability analysis and genetic correlation of traits
CPpheno = phenotype_ordered %>% mutate_at(vars('LINE'), list(as.factor))
pheid = colnames(CPpheno)[-c(1:2,grep("SampleWeight", colnames(CPpheno)))]
heritability = c()
varComp = c()

for (s in unique(CPpheno$SITE)){
  df1 = CPpheno %>% dplyr::filter(SITE==s) 
  for (phe in pheid){
    print(c(s, phe))
    if (sum(df1[,phe],na.rm = T)==0) {h2 = data.frame(SITE=s,TRAIT=phe,h2=NA,h2_SE=NA)
    Vave = cbind(data.frame(SITE=rep(s,2),TRAIT=rep(phe,2),VARIANCE=c('Va','Ve'), 
                            VarComp=rep(NA,2),VarCompSE=rep(NA,2),Zratio=rep(NA,2)))}else{
                              form =as.formula(paste(phe," ~ 1"))
                              A_phe = mmer(form, random=~vs(LINE, Gu=A), rcov=~vs(us(SITE),units),data=df1, verbose = F)
                              VaVe = cbind(data.frame(SITE=rep(s,2), TRAIT=rep(phe,2),VARIANCE=c('Va','Ve')),data.frame(summary(A_phe)$varcomp))
                              h2 = data.frame(c(s, phe,round(pin(A_phe, h2 ~ V1 / ( V1 + V2 )),2)))
                              colnames(h2) = c('SITE','TRAIT','h2','h2_SE')}
    heritability = rbind(heritability,h2)
    varComp = rbind(varComp,VaVe)}
}
###plot heritability
sites = c("PKLE",  "CLMB", "KBSM")
heritability$SITE = factor(heritability$SITE, levels = sites)

pdf('c:/Users/Li Zhang/Desktop/fourway-ionomics/1-PhenotypicAnalysis&GenstatFormatting/Heritability.pdf', onefile = T)
ggplot(heritability, aes(SITE,h2)) + geom_col(alpha=0.6) + facet_wrap(TRAIT~., nrow = 4)+theme_bw()+ ylab('Heritability')+xlab('')+
  geom_errorbar(aes(ymin=h2-h2_SE, ymax=h2+h2_SE), width=0.3)+
  theme(text = element_text(face = "bold", size = 12))+
  theme(axis.text = element_text(face = "bold", size = 8),axis.title = element_text(face = "bold", size = 10))+
  theme(strip.text = element_text(face="bold", size=12))

###plot variance components
vc = varComp
vc$SITE = factor(vc$SITE,levels=sites)
vc_va = subset(vc, VARIANCE=='Va')
vc_ve = subset(vc, VARIANCE=='Ve')
vc_ve = vc_ve %>% dplyr::mutate(VarComp=-VarComp)

ggplot(vc_va,aes(SITE,VarComp,fill=VARIANCE))+geom_col() + facet_wrap(TRAIT~.,nrow = 5,scales = 'free')+
  geom_hline(yintercept=0, color = "black",lwd=1.2)+ylab('')+xlab('')+
  geom_col(data=vc_ve,aes(SITE,VarComp,fill=VARIANCE))+facet_wrap(TRAIT~.,nrow = 5, scales = 'free')+theme_bw()+
  theme(axis.text = element_text(face = "bold", size = 6),axis.title = element_text(face = "bold", size = 10))


dev.off()
# ######
# CPpheno2 = CPpheno %>% gather('TRAITS','VALUES', SampleWeight:Cd111) %>% 
#   spread(SITE, VALUES) %>% filter(!str_detect(TRAITS, 'SampleWeight')) 
# 
# pdf('c:/Users/Li Zhang/Desktop/fourway-ionomics/1-PhenotypicAnalysis&GenstatFormatting/GeneticCorrelation.pdf', onefile=T)
# Rg <- c()
# for (i in unique(CPpheno2$TRAITS)){
#   df <- subset(CPpheno2,TRAITS==i)
#   print(i)
#   rg1 = mmer(cbind(KBSM, CLMB,PKLE)~1, random=~vs(LINE, Gu=A), rcov=~units,tolparinv = 1e-3, data = df)
#   rg= cov2cor(rg1$sigma$`u:LINE`)
#   p=corrplot::corrplot.mixed(rg, title=i,mar=c(0,0,1,0))
#   print(p)
#   
#   x1 = as.data.frame(cov2cor(rg1$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(TRAITS=i)
#   Rg = rbind(Rg,x1)
# }
# dev.off()
# write.csv(Rg,'GeneticCorrelation_among_sites.csv',row.names=F)
# 
# save.image('genetic_correlation_among_sites.RData')
# 



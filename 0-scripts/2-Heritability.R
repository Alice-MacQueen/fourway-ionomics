
### load packages
rm(list = ls())
library(tidyverse)
library(data.table)
library(pegas)
library(xlsx)
library(sommer)

## set working directory, should have genotyping file and phenotyping file
setwd(dir = '../0-data/')
getwd()

phenotype= read.csv('Ionome_cleaned.csv')
###load files needed for running sommer
load('DataFilesNeededForSommer.Rdata')

########1. heritability analysis and variance partitioning
CPpheno = phenotype %>% mutate_at(vars('LINE'), list(as.factor))
pheid = colnames(CPpheno)[-c(1:2)]
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

write.csv(heritability,'../2-Heritability&GeneticCorrelation/Heritability.csv',row.names =F )
write.csv(varComp,'../2-Heritability&GeneticCorrelation/Table_S2_VariancePartitioning.csv',row.names =F )

###2.Genetic correlation among sites for each trait
CPpheno2 = CPpheno %>% gather('TRAITS','VALUES', B:Cd) %>%spread(SITE, VALUES) 

Rg <- c()
for (i in unique(CPpheno2$TRAITS)){
  df <- subset(CPpheno2,TRAITS==i)
  print(i)
  rg1 = mmer(cbind(TX, MO,MI)~1, random=~vs(LINE, Gu=A), rcov=~units,tolparinv = 1e-3, data = df)
  rg= cov2cor(rg1$sigma$`u:LINE`)
  # p=corrplot::corrplot.mixed(rg, title=i,mar=c(0,0,1,0))
  # print(p)

  x1 = as.data.frame(cov2cor(rg1$sigma$`u:LINE`)) %>% rownames_to_column('SITE') %>% dplyr::mutate(TRAITS=i)
  Rg = rbind(Rg,x1)
}
Rg2 = Rg %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
  mutate_if(is.numeric, list(~na_if(., -Inf)))
write.csv(Rg2,'../2-Heritability&GeneticCorrelation/Table_S3_GeneticCorrelation.csv',row.names=F)


##3. test genotype by environment interaction on the trait level
lrt = NULL
for (phe in pheid){
    df1 = CPpheno
    print(phe)
    form =as.formula(paste(phe," ~ SITE"))
    ansMain = mmer(form, random=~vs(LINE, Gu=A), rcov=~units,data=df1, verbose = F)
    summary(ansMain)
    ansUS = mmer(form, random=~vs(us(SITE),LINE, Gu=A), rcov=~units,
                 data=df1, verbose = F,tolparinv = 1e-02)
    #summary(ansUS)
    x=anova(ansMain, ansUS)
    p= as.character(x$PrChisq[2])
    ps = c(phe,p)
    lrt=rbind(lrt,ps)
}
write.csv(lrt,'Likelihood-ratioTestResults.csv',row.names = F)


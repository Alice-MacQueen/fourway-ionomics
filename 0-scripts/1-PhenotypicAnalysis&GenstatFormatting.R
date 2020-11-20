
### load packages
rm(list = ls())
library(tidyverse)
library(data.table)
library(pegas)
library(xlsx)
library(rstatix)
library(ggpubr)
library(corrplot)

###0. set working directory, should have genotyping file and phenotyping file
setwd("../0-data")

###1. read the genotype .qua file to get the genotype ID
genotype = read.loci('4way.qua',skip = 4)
genotype = genotype %>% dplyr::select(id) %>% dplyr::rename(LINE=id)

###2. read the raw phenotype file and Plant ID
phenotype_raw = read.csv('Ionome.csv')
plantlist = readxl::read_excel('NSF_4WCR_Master Plant List_Final.xlsx')

###3: get the phenotypes for the grandparents
parents = phenotype_raw %>% left_join(plantlist, by=c('sample'='PLOT_GL')) %>%
  dplyr::select(SITE, LINE,B:Cd) %>% 
  filter(str_detect(LINE, "AP13|DAC|WBC|VS"))%>%
  mutate(SITE=str_replace(SITE, 'PKLE','TX'))%>%
  mutate(SITE=str_replace(SITE, 'CLMB','MO'))%>%
  mutate(SITE=str_replace(SITE, 'KBSM','MI'))%>%
  mutate_all(., list(~ifelse(.<0, NA, .)))

###4. mean and se for the grandparents for each element, and welch one-way significance test
parents_comparison = NULL
for (i in 3:20){
  print(i)
  tmp = parents[,c(1,2,i)]
  tmp2 = tmp %>% group_by(SITE, LINE) %>% 
    get_summary_stats(colnames(parents)[i], type = "mean_se")
  tmp3 = tmp2 %>% unite(X, c('mean','se'), sep = "±") %>% select(-n)%>% spread(LINE,X)
  
  tmp6 = NULL
  for (s in unique(tmp3$SITE)){
    if (i==15|i==16){tmp4 = tmp4 %>% filter(LINE!='VS16')}else{tmp4 = tmp %>% filter(SITE==s)}
    colnames(tmp4)[3] = "ion"
    res.aov <- tmp4 %>% anova_test(ion ~ LINE)###anova test, assuming equal variance
    res.oneway = oneway.test(tmp4$ion~tmp4$LINE, na.action = na.omit ) ##one way test,not necessarity equal variance
    p_aov = res.aov$p
    p_one = res.oneway$p.value 
    tmp5 = c(s, p_aov, p_one)
    tmp6 = rbind(tmp6, tmp5)
    colnames(tmp6)=c('SITE',"P-aov","P-oneway")
  }
  tmp7 = tmp3 %>% left_join(tmp6 %>% as.data.frame())
  parents_comparison = rbind(parents_comparison,tmp7)
}

write.csv(parents_comparison,'../1-PhenotypicAnalysis&GenstatFormatting/Table1_ParentsComparisonForEachIon.csv',row.names = F)

parents_comparison2 = parents_comparison %>% 
  separate(AP13,c('AP13_mean','AP13_se'),sep="±")%>%
  separate(DAC,c('DAC_mean','DAC_se'),sep="±")%>%
  separate(VS16,c('VS16_mean','VS16_se'),sep="±")%>%
  separate(WBC,c('WBC_mean','WBC_se'),sep="±")
write.csv(parents_comparison2,'../1-PhenotypicAnalysis&GenstatFormatting/Table1_ParentsComparisonForEachIon2.csv',row.names = F)

##5. formating to genstat and get the mean value for each line
phenotype_avg = phenotype_raw %>% left_join(plantlist, by=c('sample'='PLOT_GL')) %>%
  dplyr::select(SITE, LINE,B:Cd) %>% group_by(SITE, LINE) %>%
  summarise_all(list(mean),na.rm=T) %>% ungroup()

phenotype_ordered = vector('list', 0)
for (s in unique(phenotype_avg$SITE)){
  phenotype_ordered[[s]] = phenotype_avg %>%  filter(SITE==s) %>% right_join(genotype)%>% mutate(SITE = s)
}
phenotype_ordered =bind_rows(phenotype_ordered)

###removing negative values is traits-based, 
phenotype_positive = phenotype_ordered %>%
  mutate_all(., list(~ifelse(.<0, NA, .))) %>% 
  gather(traits, values,-SITE, -LINE)

phenotype_clean = NULL
for (phe in unique(phenotype_positive$traits)){
  for (s in unique(phenotype_positive$SITE)){
  tmp = phenotype_positive %>% filter(traits==phe)%>% filter(SITE==s)
  outliers = boxplot.stats(tmp$values,coef = 2)$out
  tmp$values[tmp$values %in% outliers] = NA
  LineRemove = tmp %>% group_by(LINE) %>% summarise(Count=n(),SUM = sum(values, na.rm=T)) %>% filter(SUM==0) 
  tmp2  = tmp %>% filter(!LINE %in% LineRemove$LINE)
  phenotype_clean = rbind(phenotype_clean,tmp2)
  }
}
phenotype_clean = phenotype_clean %>% mutate(SITE=str_replace(SITE, 'PKLE','TX'))%>%
  mutate(SITE=str_replace(SITE, 'CLMB','MO'))%>%
  mutate(SITE=str_replace(SITE, 'KBSM','MI'))


###6: check if ion differ between sites, welch one-way test
tmp5 = NULL
tmp6 = NULL
for (i in unique(phenotype_clean$traits)){
  tmp = phenotype_clean %>% filter(traits==i)
  colnames(tmp)[4]=i
  tmp2 = tmp %>% group_by(SITE) %>% 
    get_summary_stats(colnames(tmp)[4], type = "mean_se")
  tmp3 = tmp2 %>% unite(X, c('mean','se'), sep = "±") %>% select(-n)%>% spread(SITE,X)
  
    colnames(tmp)[4] = "ion"
    res.aov <- tmp %>% anova_test(ion ~ SITE)###anova test, assuming equal variance
    res.oneway = oneway.test(tmp$ion~tmp$SITE) ##one way test,not necessarity equal variance
    p_aov = res.aov$p
    p_one = res.oneway$p.value 
    tmp4 = c(i, p_aov, p_one)
    tmp5 = rbind( tmp5,tmp4)
    colnames(tmp5)=c('variable',"P-aov","P-oneway")
    tmp6 = rbind(tmp6, tmp3)
  }
  site_comparison = tmp6 %>% left_join(tmp5%>% as.data.frame())

  write.csv(site_comparison,'../1-PhenotypicAnalysis&GenstatFormatting/Table2_SiteComparisonForEachIon.csv',row.names = F)
  
 site_comparison2 = site_comparison %>% 
    separate(MI,c('MI_mean','MI_se'),sep="±")%>%
    separate(MO,c('MO_mean','MO_se'),sep="±")%>%
    separate(TX,c('TX_mean','TX_se'),sep="±")
  write.csv(site_comparison2,'../1-PhenotypicAnalysis&GenstatFormatting/Table2_SitesComparisonForEachIon2.csv',row.names = F)
  
# ###7: check if ion differs between pairs of sites using ks.test, not used
#   KS_Pvalues = NULL
#   for (i in unique(phenotype_clean$traits)){
#     tmp = phenotype_clean %>% filter(traits==i)%>% spread(SITE,values)
#     MO_MI = ks.test(tmp$MO, tmp$MI)
#     MO_TX = ks.test(tmp$MO, tmp$TX)
#     MI_TX = ks.test(tmp$MI, tmp$TX)
#     MO_MI_P = MO_MI$p.value
#     MO_TX_P = MO_TX$p.value
#     MI_TX_P = MI_TX$p.value
#     ks_p = cbind(i, MO_MI_P, MO_TX_P,MI_TX_P)
#     KS_Pvalues = rbind(KS_Pvalues, ks_p)
#   }
# write.csv(KS_Pvalues,'C:/Users/Li Zhang/Desktop/fourway-ionomics/1-PhenotypicAnalysis&GenstatFormatting/Table2_Pvalues_ks_test.csv',row.names = T)
#   

phenotype_clean2 = phenotype_clean %>% spread(traits, values)
write.csv(phenotype_clean2,'Ionome_cleaned.csv',row.names = F)

###8:phenotypic correlation,and rearrange elements a little bit
phenotype_clean2 = phenotype_clean2 %>% dplyr::select(SITE, LINE, P, K, Ca, Mg,     ##macronutrients
                                                      Rb, Sr,                       ## analogues
                                                      B, Mn, Fe, Co, Cu, Zn, Se, Mo,## micronutrients
                                                      Na, Al, As, Cd)               ## harmful

Correlations = NULL
for(s in unique(phenotype_clean2$SITE)){
  tmp = phenotype_clean2 %>% filter(SITE==s)
  mcor = round(cor(tmp[,-c(1,2)],use = "complete.obs"),2)
  ##hide upper triangle
  upper = mcor
  upper[upper.tri(mcor)] = ""
  upper = as.data.frame(upper)
  corr = upper  %>% rownames_to_column('Element') %>% mutate('SITE'=s)
  Correlations = rbind(Correlations, corr)
}
write.csv(Correlations, '../1-PhenotypicAnalysis&GenstatFormatting/Table_S1_PhenotypicCorrelationsForEachElementAtEachSite.csv',row.names = F)

###9: formatting for Genstat
#######number of genotypes for each traits varies after outlier removal, so need to build genstat file for each trait.
setwd(dir = paste(dir_main, '2-GenstatRunning',sep = '/'))
getwd()
geno=read.loci('4way_reduced_marker.loc',header = F)
colnames(geno) = c('Marker','Phase','NoIdea',as.character(genotype$LINE))

for (phe in unique(phenotype_positive$traits)){
  tmp = phenotype_positive %>% filter(traits==phe)
  outliers = boxplot.stats(tmp$values,coef = 2)$out
  #outliers
  tmp$values[tmp$values %in% outliers] = NA
  LineRemove = tmp %>% group_by(LINE) %>% summarise(Count=n(),SUM = sum(values, na.rm=T)) %>% filter(SUM==0) 
  print(length(LineRemove$LINE))
  Genstat = tmp %>% filter(!LINE %in% LineRemove$LINE) %>% mutate(SITE2 = SITE)
  colnames(Genstat)[4] = phe
  Genstat[is.na(Genstat)] <- "*"
  write.csv(Genstat, paste0('GenstatFormated_',phe,'.csv'),row.names = F)
  
  ###read the loci file and remove the genotype information for these lines which had no phenotyping data

  geno_sub = geno[,!colnames(geno) %in% LineRemove$LINE]
  write.loci(geno_sub,'geno_subset.loc',loci.sep = '\t',quote=F,row.name=F,col.names=F)
  
  Genotype_sub= data.frame(LINE=as.numeric(colnames(geno_sub)[-c(1:3)]))
  write.csv(Genotype_sub,'Genotype_sub.csv',row.names = F)
  
  ###6. put the loci and genotype_sub file together for Genstat format
  y = readLines('geno_subset.loc')
  z = readLines('Genotype_sub.csv')
  
  cat('name = cross\npopt = CP\nnloc=738',file=paste0('4way_sub_',phe,'.loc'),sep = '\n',append = F)
  cat(paste0('nind = ',length(Genotype_sub$LINE)),file=paste0('4way_sub_',phe,'.loc'),sep = '\n',append = T)
  cat(y,file=paste0('4way_sub_',phe,'.loc'),sep = '\n',append = T)
  cat('individual names:',file=paste0('4way_sub_',phe,'.loc'),sep = '\n',append = T)
  cat(z[-1],file=paste0('4way_sub_',phe,'.loc'),sep = '\n',append = T)
  ###END
  file.remove('geno_subset.loc','Genotype_sub.csv')
  ### use '4way_sub_phe.loc', 'GenstatFormated_phe.csv', '4way_male_reduced_marker.map' file in Genstat
  ### then run Genstat
}





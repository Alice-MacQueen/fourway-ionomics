
### load packages
rm(list = ls())
library(tidyverse)
library(data.table)
library(pegas)
library(xlsx)
library(corrplot)

###0. set working directory, should have genotyping file and phenotyping file
dir_main = 'c://Users/Li Zhang/Desktop/IonomicsData/'
setwd(dir = paste(dir_main,'0-data',sep = '/'))
getwd()

###1. read the genotype .qua file to get the genotype ID
genotype = read.loci('4way.qua',skip = 4)
genotype = genotype %>% dplyr::select(id) %>% dplyr::rename(LINE=id)

###2. read the raw phenotype file
phenotype_raw = read.csv('Ionome.csv')
plantlist = readxl::read_excel('NSF_4WCR_Master Plant List_Final.xlsx')

##3. formating to genstat and get the mean value for each line
phenotype_avg = phenotype_raw %>% left_join(plantlist, by=c('sample'='PLOT_GL')) %>%
  dplyr::select(SITE, LINE,SampleWeight:Cd111) %>% group_by(SITE, LINE) %>%
  summarise_all(list(mean),na.rm=T) %>% ungroup()

phenotype_ordered = vector('list', 0)
for (s in unique(phenotype_avg$SITE)){
  phenotype_ordered[[s]] = phenotype_avg %>%  filter(SITE==s) %>% right_join(genotype)%>% mutate(SITE = s)
}
phenotype_ordered =bind_rows(phenotype_ordered)
###plot the historgram for each trait, check for negative value and outliers
out1 = phenotype_ordered %>% gather(traits, values,-SITE, -LINE)
pdf(paste(dir_main,'1-PhenotypicAnalysis&GenstatFormatting/Histograms of Iron Traits_before removing outliers.pdf', sep = '/'),onefile = T)
par(mfrow=c(2,2))
for (phe in unique(out1$traits)){
  x = out1 %>% filter(traits==phe)
  hist(x$values, main = phe, xlab='')
}
dev.off()

###removing negative values is traits-based, some needs to be set as 0 from negative values, some needs to be set as NA
phenotype_positive = phenotype_ordered %>%
  mutate_all(., list(~ifelse(.<0, NA, .))) %>% 
  gather(traits, values,-SITE, -LINE)

pdf(paste(dir_main,'1-PhenotypicAnalysis&GenstatFormatting/Histograms of Iron Traits_after removing outliers.pdf',sep='/'),onefile = T)
par(mfrow=c(2,2))
for (phe in unique(phenotype_positive$traits)){
  tmp = phenotype_positive %>% filter(traits==phe)
  #hist(tmp$values, main = phe, xlab='')
  #summary(tmp$values)
  outliers = boxplot.stats(tmp$values,coef = 2)$out
  #outliers
  tmp$values[tmp$values %in% outliers] = NA
  #summary(tmp$values)
  hist(tmp$values, main = phe, xlab='')
}
dev.off()

###correlation plots among traits
phenotype_clean = NULL
for (phe in unique(phenotype_positive$traits)){
  tmp = phenotype_positive %>% filter(traits==phe)
  LineRemove = tmp %>% group_by(LINE) %>% summarise(Count=n(),SUM = sum(values, na.rm=T)) %>% filter(SUM==0) 
  tmp2  = tmp %>% filter(!LINE %in% LineRemove$LINE)
  phenotype_clean = rbind(phenotype_clean,tmp2)
  }
phenotype_clean = phenotype_clean %>% spread(traits, values)

pdf(paste0(dir_main,'1-PhenotypicAnalysis&GenstatFormatting/PhenotypicCorrelation.pdf'), onefile = T, width = 10, height=10)
M = cor(phenotype_clean[,-c(1:2)], use = 'pairwise.complete.obs')
corrplot::corrplot.mixed(M,lower.col = 'blue')
mtext('Correlation Across Sites',side=2)

for(s in unique(phenotype_clean$SITE)){
  tmp = phenotype_clean %>% filter(SITE==s)
  M = cor(tmp[,-c(1:2)], use = 'pairwise.complete.obs')
  corrplot::corrplot.mixed(M,lower.col = 'blue')
  mtext(paste0('Correlation at_',s), side=2)
}
dev.off()

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





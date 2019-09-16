### Formatting Genstat Results##
rm(list=ls())
library(tidyverse)
library(reshape2)
library(xlsx)
library(ape)
setwd("c://Users/Li Zhang/Desktop/IonomicsData/4-CandidateGeneSearch/")
### Read QTL info
QTL_all = read.csv('c:/Users/Li Zhang/Desktop/IonomicsData/3-GenstatOutput/QTLlistForAllTraits_flankmarker.csv')

QTL = QTL_all %>% dplyr::mutate_at(vars(down_marker, up_marker, TRAIT), list(as.character)) %>% 
  separate(down_marker, c('chr','pos_lo'),sep = "_") %>% separate(up_marker, c('chr_hi','pos_hi'), sep="_") %>% 
  dplyr::mutate_at(vars(pos_lo, pos_hi), list(as.numeric)) %>% rowwise()%>%
  mutate(pos_lo= pos_lo*1000000, pos_hi=pos_hi*1000000)%>%
  dplyr:: select(chr,QTL, pos_lo, pos_hi, TRAIT) 


##Read annotation file
genesid = read.gff('Pvirgatum_516_v5.1.gene.gff3.gz')
genesid2 = genesid %>% filter(type=='gene') %>% dplyr::rename(ID=attributes) %>%  separate(ID, c('A','B','C'), sep="=|\\;|\\.") %>%
  unite('locusName','B','C',sep = '.') %>%  dplyr::select(seqid, source, type, start, end, locusName) 

annot = read.delim2('Pvirgatum_516_v5.1.annotation_info.txt',sep = "\t",header = T)
annot1 = annot %>% dplyr::select(locusName, Best.hit.arabi.name:Best.hit.rice.name, rice.defline) %>% unique()

###the annotation file with gene ids
annot_gene = genesid2 %>% left_join(annot1)

Counts = NULL
CandidateGeneAll = NULL
for (i in 1:nrow(QTL)){ 
  gene_id_ind = NULL
  print(i)
  
  annot_gene_sub = annot_gene %>% dplyr::filter(seqid==QTL$chr[i])
  for (j in 1:nrow(annot_gene_sub)){
    # print(c(i,j))
    if ((annot_gene_sub[j,4] >= as.numeric(QTL$pos_lo[i])) &  (annot_gene_sub[j,5]<= as.numeric(QTL$pos_hi[i])))
    {gene_id_ind = rbind(gene_id_ind,annot_gene_sub[j,])}
  }
  gene_id_ind2 = gene_id_ind %>% dplyr::mutate_at(vars(seqid),list(as.character)) %>% dplyr::rename(chr=seqid) %>% 
    dplyr::select(locusName,start, end, "Best.hit.arabi.name","arabi.symbol","arabi.defline","Best.hit.rice.name" , "rice.defline" ) %>%
    separate(Best.hit.arabi.name,c('Best.hit.arabi.name','ver'),sep = '\\.') %>% separate(Best.hit.rice.name, c('Best.hit.rice.name','ver2'),sep = '\\.') %>%
    dplyr::select(-ver, -ver2)%>% group_by(locusName) %>% slice(1) %>% ungroup()
  
  QTL_sub=bind_rows(replicate(nrow(gene_id_ind2), QTL[i,], simplify = FALSE))
  gene_id_ind2 = bind_cols(QTL_sub, gene_id_ind2)
  
  counts = data.frame(QTL=QTL$QTL[i], pvGenes=length(gene_id_ind2$locusName), ArabGenes=length(unique(gene_id_ind2$Best.hit.arabi.name))-1,
                      RiceGenes=length(unique(gene_id_ind2$Best.hit.rice.name))-1)
  Counts = rbind(Counts, counts)
  
  write.csv(gene_id_ind2,paste0('CandidateGeneLists_',QTL$QTL[i],'_',QTL$TRAIT[i], '.csv'), row.names = F)
  CandidateGeneAll = rbind(CandidateGeneAll,gene_id_ind2)
}
write.csv(Counts,'NumberOfCandidateGenesSummary.csv',row.names = F)

write.csv(CandidateGeneAll,'CandidateGene_ALl.csv',row.names = F)


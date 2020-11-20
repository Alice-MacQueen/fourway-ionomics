### Formatting Genstat Results##
rm(list=ls())
library(tidyverse)
library(reshape2)
library(xlsx)
library(ape)
options(warn=-1)
setwd("c://Users/Li Zhang/Desktop/fourway-ionomics/4-CandidateGeneSearch/")
### Read QTL info
QTL_all = read.csv('c:/Users/Li Zhang/Desktop/fourway-ionomics/3-GenstatOutput/QTLlistForAllTraits_flankmarker.csv')

QTL = QTL_all %>% dplyr::mutate_at(vars(down_marker, up_marker, TRAIT), list(as.character)) %>% 
  separate(down_marker, c('chr','pos_lo'),sep = "_") %>% separate(up_marker, c('chr_hi','pos_hi'), sep="_") %>% 
  dplyr::mutate_at(vars(pos_lo, pos_hi), list(as.numeric)) %>% rowwise()%>%
  mutate(pos_lo= pos_lo*1000000, pos_hi=pos_hi*1000000)%>%
  dplyr:: select(chr,QTL, pos_lo, pos_hi, TRAIT) 

RiceGene = read.csv('c:/Users/Li Zhang/Desktop/fourway-ionomics/0-data/Fun.Rice.Genes.3427.csv')

####Read gff3 file and annotation file
gff3 = read.gff('c:/Users/Li Zhang/Desktop/fourway-ionomics/0-data/Pvirgatum_516_v5.1.gene.gff3.gz')
genes = gff3 %>% filter(type=='gene') %>% dplyr::rename(ID=attributes) %>%  separate(ID, c('A','B','C'), sep="=|\\;|\\.") %>%
  unite('locusName','B','C',sep = '.') %>%  dplyr::select(seqid, start, end, locusName) 

annot = read.delim2('c:/Users/Li Zhang/Desktop/fourway-ionomics/0-data/Pvirgatum_516_v5.1.annotation_info.txt',sep = "\t",header = T)
##wroking on the annotation file 
annot2 = annot %>% dplyr::select(locusName, GO, "Best.hit.arabi.name","arabi.symbol","arabi.defline","Best.hit.rice.name" , "rice.defline" )%>%
  separate(Best.hit.arabi.name,c('Best.hit.arabi.name','ver'),sep = '\\.') %>% separate(Best.hit.rice.name, c('Best.hit.rice.name','ver2'),sep = '\\.') %>%
  dplyr::select(-ver, -ver2) %>% distinct(locusName,.keep_all = T)

###write the annotation file for enrichment analysis
geneID2GO = annot2 %>% dplyr::select(locusName, GO) %>% dplyr::filter(GO != "") 
write.table(geneID2GO, file ='annotation_for_enrichment.txt',sep=" ", quote = F, row.names = F, col.names = F )
####

genes_annot = genes %>% left_join(annot2) %>% left_join(RiceGene, by=c('Best.hit.rice.name'='OsID'))


#######get the candidate genes for each QTL intervel for the traits
Counts = NULL
CandidateGeneAll = NULL
for (i in 1:nrow(QTL)){ 
  gene_list = genes_annot%>% dplyr::filter(seqid==QTL$chr[i]) %>% 
    filter(start >= as.numeric(QTL$pos_lo[i]) & end <= as.numeric(QTL$pos_hi[i]))
  
  QTL_sub=bind_rows(replicate(nrow(gene_list), QTL[i,], simplify = FALSE))
  gene_list = bind_cols(QTL_sub, gene_list)
  
  counts = data.frame(QTL=QTL$QTL[i], TRAIT=QTL$TRAIT[i], pvGenes=length(gene_list$locusName), ArabGenes=length(unique(gene_list$Best.hit.arabi.name))-1,
                      RiceGenes=length(unique(gene_list$Best.hit.rice.name))-1)
  Counts = rbind(Counts, counts)
  #write.csv(gene_list,paste0('CandidateGeneLists_',QTL$QTL[i],'_',QTL$TRAIT[i], '.csv'), row.names = F)
  CandidateGeneAll = rbind(CandidateGeneAll,gene_list)
}
write.csv(Counts,'NumberOfCandidateGenesSummary.csv',row.names = F)
write.csv(CandidateGeneAll,'CandidateGene_ALL.csv',row.names = F)


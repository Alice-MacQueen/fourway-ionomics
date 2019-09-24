### Formatting Genstat Results##
rm(list=ls())
library(tidyverse)
library(reshape2)
library(xlsx)
library(ape)

setwd("c://Users/Li Zhang/Desktop/fourway-ionomics/5-EnrichmentAnalysis/")

##Read gff3 file and annotation file
gff3 = read.gff('c:/Users/Li Zhang/Desktop/fourway-ionomics/0-data/Pvirgatum_516_v5.1.gene.gff3.gz')
genes = gff3 %>% filter(type=='gene') %>% dplyr::rename(ID=attributes) %>%  separate(ID, c('A','B','C'), sep="=|\\;|\\.") %>%
  unite('locusName','B','C',sep = '.') %>%  dplyr::select(seqid, source, type, start, end, locusName) 

annot = read.delim2('c:/Users/Li Zhang/Desktop/fourway-ionomics/0-data/Pvirgatum_516_v5.1.annotation_info.txt',sep = "\t",header = T)

genes_annot = genes %>% left_join(annot) %>% distinct(locusName,.keep_all = T)%>% dplyr::select(seqid, start, end, locusName, GO, "Best.hit.arabi.name","arabi.symbol","arabi.defline","Best.hit.rice.name" , "rice.defline" )%>%
  separate(Best.hit.arabi.name,c('Best.hit.arabi.name','ver'),sep = '\\.') %>% separate(Best.hit.rice.name, c('Best.hit.rice.name','ver2'),sep = '\\.') %>%
  dplyr::select(-ver, -ver2)%>% group_by(locusName) %>% slice(1) %>% ungroup()

##read candidate gene list for all traits
CandidateGeneAll = read.csv('c:/Users/Li Zhang/Desktop/fourway-ionomics/4-CandidateGeneSearch/CandidateGene_ALl.csv')

gene_in_genome = genes_annot %>% filter(GO !="")

go_in_genome = gene_in_genome %>% dplyr::select(locusName, GO) %>% 
  separate(GO, into = c("GO1", "GO2", "GO3", "GO4", "GO5", "GO6", "GO7", "GO8","GO9", "GO10"), sep = ",") %>%
  gather('GOTerm','GO', -locusName) %>% dplyr::select(-GOTerm)%>% filter(!is.na(GO))%>%  group_by(GO) %>% summarise(count = n()) %>%
  arrange(desc(count)) 

###go in samples
gene_in_sample = CandidateGeneAll %>% filter(GO !="") %>% group_by( TRAIT,locusName) %>%summarise(count = n()) %>% ungroup()

go_in_sample= CandidateGeneAll %>% filter(GO !="") %>% dplyr::select(locusName,TRAIT, GO) %>%
  separate(GO, into = c("GO1", "GO2", "GO3", "GO4", "GO5", "GO6", "GO7", "GO8","GO9", "GO10"), sep = ",") %>%
  gather('GOTerm','GO', -locusName, -TRAIT) %>% dplyr::select(-GOTerm)%>% group_by(TRAIT,GO) %>% summarise(count = n())%>%
  arrange(TRAIT, desc(count))%>% filter(str_detect(GO,'GO'))


go_significant = NULL
###find the number of candidate genes for each trait
for (tr in unique(CandidateGeneAll$TRAIT)){
  gene_in_sample_trait = gene_in_sample %>% filter(TRAIT==tr)
  go_in_sample_trait = go_in_sample %>% filter(TRAIT==tr) 
  
  for (i in 1:nrow(go_in_sample_trait)){
    go_sample = as.numeric(go_in_sample_trait[i,"count"])
    go_genome = go_in_genome %>% filter(GO==go_in_sample_trait$GO[i]) %>% dplyr::select(count)%>% as.numeric()
    
    
    ###build the contigency table
    GOmatrix = matrix(c(go_sample,(go_genome - go_sample), 
                        (nrow(gene_in_sample_trait)- go_sample), 
                        (nrow(gene_in_genome)-nrow(gene_in_sample_trait)- (go_genome - go_sample))), 
                      nrow = 2, 
                      dimnames = list(c("GO in Trait", "Go not in Trait"), 
                                      c("Gene in GO", "Genes not in GO")))
    
    go_in_sample_trait[i, "p.value"] = fisher.test(GOmatrix)$p.value
  }
  cutoff = quantile(go_in_sample_trait$p.value, probs = 0.05)
  go_significant_trait = go_in_sample_trait %>% filter (p.value<cutoff)
  go_significant = rbind(go_significant, go_significant_trait)
  print(c(tr,i,cutoff, nrow(go_in_sample_trait)))
}
write.csv(go_significant,'SignificantGOtermsForTraits.csv',row.names = F)
length(unique(go_significant$GO))

########get the function for these GO terms
library(ontologyIndex)
ontology = get_OBO('c:/Users/Li Zhang/Desktop/fourway-ionomics/0-data/go-basic.obo', extract_tags = 'everything')
terms = unique(go_significant$GO)
go_id = data.frame(GO=ontology$id[terms])
go_name = data.frame(GO=ontology$id[terms], NAME=ontology$name[terms])
go_namespace = data.frame(FUNCTION=unlist(ontology$namespace[terms]))
go_namespace = go_namespace %>% rownames_to_column('GO')

go_functions = go_id %>% left_join(go_name) %>% left_join(go_namespace) %>% filter(complete.cases(.))
write.csv(go_functions, 'FunctionsOfSignificantGOterms.csv',row.names = F)

go_functions_trait = go_significant %>% left_join(go_functions) %>% dplyr::select(-count, -'p.value')
write.csv(go_functions_trait,'FunctionsOfSignificantGOtermsForEachTrait.csv',row.names = F)

go_genes = CandidateGeneAll %>% filter(GO !="") %>% 
  dplyr::select(-chr, -pos_lo, -pos_hi, -seqid, -start,-end) %>%
  separate(GO, into = c("GO1", "GO2", "GO3", "GO4", "GO5", "GO6", "GO7", "GO8","GO9", "GO10"), sep = ",") %>%
  gather('GOTerm','GO', GO1:GO10) %>% dplyr::select(-GOTerm)

go_genes =go_functions_trait %>% left_join(go_genes)
write.csv(go_genes,'Functions&CandidateGenesForSignificantGOterms.csv',row.names = F)



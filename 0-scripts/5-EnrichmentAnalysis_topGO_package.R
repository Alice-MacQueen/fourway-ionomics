### Formatting Genstat Results##
rm(list=ls())
library(tidyverse)
library(topGO)

setwd("c://Users/Li Zhang/Desktop/fourway-ionomics/0-data/")
##read the annotation file that specifically for enrichemnt analysis
geneID2GO = readMappings(file = "annotation_for_enrichment.txt", sep = " ")  

geneUniverse = names(geneID2GO)

CandidateGeneAll = read.csv('c:/Users/Li Zhang/Desktop/fourway-ionomics/4-CandidateGeneSearch/CandidateGene_ALl.csv')

CandidateGeneAll = CandidateGeneAll %>% dplyr::select(TRAIT,locusName, GO) %>% dplyr::filter(GO!="")

go_significant_all_trait = NULL
for (tr in unique(CandidateGeneAll$TRAIT)){
   tmp = CandidateGeneAll %>% filter(TRAIT==tr) %>% dplyr::select(locusName)
   genesOfInterest <- as.character(tmp$locusName) 
   geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
   names(geneList) <- geneUniverse
   head(geneList)
   table(geneList)

   ###test for BP, MF, CC:
   go_significant = NULL
   for (pr in c("BP","MF","CC")){
       GOdata = new("topGOdata",
             ontology = pr,
             allGenes = geneList,
             annot = annFUN.gene2GO, gene2GO = geneID2GO)
       resultFisher= runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  
       tab = GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
               numChar = 120)
       tab_significant = tab %>% filter (raw.p.value < 0.01) %>% mutate(Trait=tr, pr)
       go_significant = rbind(tab_significant, go_significant)
   }
   go_significant_all_trait = rbind(go_significant_all_trait, go_significant)
}

write.csv(go_significant_all_trait, "../5-EnrichmentAnalysis/SignificantGOsForAllTraits.csv", row.names = F)

BiocManager::install('grimbough/biomaRt')
library(biomaRt)

#function to retrieve correspondence between zebra finch genes and human Ensembl IDs (barn swallow was not on the website at the time fo the analysis)
getOrthos <- function(input_org,output_org,one2one) { 
     tmp_input <- paste(input_org,"_gene_ensembl",sep = "")
     tmp_mart <- useMart("ensembl", dataset = tmp_input,host = "www.ensembl.org")
     tmp_attr <- listAttributes(tmp_mart)
     tmp_genes <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"), mart = tmp_mart)
     tmp_genes <- tmp_genes[grep("protein_coding",tmp_genes$gene_biotype),]
     attr_names <- tmp_attr[grep(output_org,tmp_attr$name),"name"]
     tmp_orthos <- getBM(attributes = c("ensembl_gene_id","external_gene_name",attr_names), filters = "ensembl_gene_id",values = tmp_genes$ensembl_gene_id, mart = tmp_mart)  
     
     if(one2one == TRUE){
     tmp_orthos <- tmp_orthos[grep("ortholog_one2one",tmp_orthos[,paste(output_org,"_homolog_orthology_type",sep = "")]),]
     }  
  if(one2one == FALSE){
      tmp_orthos <- tmp_orthos
       }    
    tmp_orthos
  }

convert_table <- getOrthos("tguttata","hsapiens",TRUE) #makes table converting bird gene_id to hsap_enembl_id.
write.table(convert_table, file="tguttata_hsampies.txt", sep = "\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
#from this generate a txt file with only the "Universe genes" and the "gage_expr*.txt" files

#table with human genes and Go terms
bm <- useMart("ensembl")                                         
bm <- useDataset("hsapiens_gene_ensembl", mart=bm)
GO_tbl <- getBM(mart=bm, attributes=c('ensembl_gene_id','go_id','name_1006')) #human go list
write.table(GO_tbl, file="GO_tbl.txt", sep = "\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
#comvert into csv

library(tidyverse)
library(gage)

#load the universe genes
Universe_genes <- read.table("Universe_genes.txt", header=FALSE) #list of human ens ids in universe (that maps to Zebra finch IDs)
#load the Human IDs and GO terms
GO_tbl <- read.csv2("GO_tbl.csv",header=TRUE) %>% as_tibble
#mantain only the genes from our universe
GO_tbl2 <- GO_tbl %>% filter(ensembl_gene_id %in% Universe_genes$V1)

#generate nested list with a list of genes for each GO term
g_sets <- list()
Unique_go <- unique(GO_tbl2$go_id)

for (go in Unique_go){
  GO_tbl_subset <- GO_tbl2 %>% filter(go_id == go)
  tmp_genes <- GO_tbl_subset$ensembl_gene_id
  g_sets[[go]] <- tmp_genes
}

#gage analysis
gage_exprs_606CONS <- read.table("gage_exprs_606CONS.txt", row.names = 1)
gage_exprs_606ACC <- read.table("gage_exprs_606ACC.txt", row.names = 1)

gage_results_606CONS <- gage(gage_exprs_606CONS,g_sets)$greater
gage_results_606ACC <- gage(gage_exprs_606ACC,g_sets)$greater

write.table(gage_results_606CONS, file = "gage_results_606CONS.txt",sep = "\t")
write.table(gage_results_606ACC, file = "gage_results_606ACC.txt",sep = "\t")











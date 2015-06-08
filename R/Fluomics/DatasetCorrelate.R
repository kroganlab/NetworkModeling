source('R/DB/DBFunctions.R')
library(reshape2)
library(pheatmap)
library(RColorBrewer)

absmax = function(x){
  if(abs(min(x)) > abs(max(x))){
    min(x)
  }else{
    max(x)  
  }
}

CELL_LINE = 'MOUSE_LUNG'
STRAIN = 'H1N1'
LFC = 2
FDR = 0.0001

###############
## GET DATASETS

CONDITION='.*'
genomics_mouse_H1N1 = getGenomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'genomics_rnaseq_mouse', condition_1 = '.*', condition_2 = '.*')
prot_ph_mouse_H1N1 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ph_sites_mouse', mapping = 'mygene_mouse_1to1', condition_1 = '.*', condition_2 = '.*')
prot_ub_mouse_H1N1 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ub_sites_mouse', mapping = 'mygene_mouse_1to1', condition_1 = '.*', condition_2 = '.*')
metabol_mouse_H1N1_known_only = getMouseMetabolomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, condition_1 = '.*', condition_2 = '.*', only_known_metabolites = T, mouse_human_orthology_confidence = 1)
mouse_H1N1 = rbind(genomics_mouse_H1N1, prot_ph_mouse_H1N1, prot_ub_mouse_H1N1, metabol_mouse_H1N1_known_only)


## get all these genes from the DB, regardless of cutoff (union of unique genes in at least one conditions)
mouse_H1N1_entrez_unique = unique(mouse_H1N1$entrez_id)
na.omit(mouse_H1N1_entrez_unique)

genomics_mouse_H1N1_consensus = getGenomics(cell_line = CELL_LINE, strain = STRAIN, table = 'genomics_rnaseq_mouse', condition_1 = '.*', condition_2 = '.*', by = "entrez_id", id_list = mouse_H1N1_entrez_unique)
prot_ph_mouse_H1N1_consensus = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, table = 'proteomics_ph_sites_mouse', mapping = 'mygene_mouse_1to1', condition_1 = '.*', condition_2 = '.*', by="entrez_id", id_list = mouse_H1N1_entrez_unique)
prot_ub_mouse_H1N1_consensus = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, table = 'proteomics_ub_sites_mouse', mapping = 'mygene_mouse_1to1', condition_1 = '.*', condition_2 = '.*', by="entrez_id", id_list = mouse_H1N1_entrez_unique)
prot_mb_mouse_H1N1_consensus = getMouseMetabolomics(cell_line = CELL_LINE, strain = STRAIN, condition_1 = '.*', condition_2 = '.*', by="entrez_id", id_list = mouse_H1N1_entrez_unique, only_known_metabolites = T)
mouse_H1N1_consensus = rbind(genomics_mouse_H1N1_consensus, prot_ph_mouse_H1N1_consensus, prot_ub_mouse_H1N1_consensus, prot_mb_mouse_H1N1_consensus)


###################
## MAKE DATA MATRIX

## convert q-values to log10 q-values
mouse_H1N1_consensus = data.frame(mouse_H1N1_consensus, log10_qvalue=-log10(mouse_H1N1_consensus$q_value))
mouse_H1N1_consensus$log10_qvalue[is.infinite(mouse_H1N1_consensus$log10_qvalue)] = 16
## data matrix with p-values
mouse_H1N1_consensus_w = dcast(entrez_id ~ condition_2+omics_type, data=mouse_H1N1_consensus, value.var='log10_qvalue', fun.aggregate = max)
## data matrix with fold-changes
mouse_H1N1_consensus_w = dcast(entrez_id ~ condition_2+omics_type, data=mouse_H1N1_consensus, value.var='log2fc', fun.aggregate = absmax)

#########################
## get gene names from DB

con = myConnect()
res = dbGetQuery(conn = con, statement = sprintf("select entrez_id_mygene, gene_name from mygene_mouse where entrez_id_mygene regexp '%s' group by entrez_id_mygene, gene_name", paste(mouse_H1N1_consensus_w$entrez_id, collapse='|')))
dbDisconnect(con)
mouse_H1N1_consensus_w = merge(mouse_H1N1_consensus_w, res[,c('entrez_id_mygene','gene_name')], by.x='entrez_id', by.y='entrez_id_mygene')
mouse_H1N1_consensus_w$gene_name = paste(mouse_H1N1_consensus_w$gene_name, mouse_H1N1_consensus_w$entrez_id, sep=' - ')

#########################
## convert for heatmap 

rownames(mouse_H1N1_consensus_w) = mouse_H1N1_consensus_w$gene_name
mouse_H1N1_consensus_w = as.matrix(mouse_H1N1_consensus_w[,-c(1,ncol(mouse_H1N1_consensus_w))])
mouse_H1N1_consensus_w[is.infinite(mouse_H1N1_consensus_w)]=0
pheatmap(mouse_H1N1_consensus_w, file='~/Desktop/mouse_H1N1.pdf', cellheight = 10, cellwidth=10, clustering_distance_cols = 'correlation', clustering_distance_rows = 'correlation', color = c(rev(brewer.pal(n=3,name = 'Blues')), '#FFFFFF','#FFFFFF',brewer.pal(n=3,name = 'Reds')), breaks = c(-10, -6, -4,-2, 0, 2, 4, 6, 10), legend_labels = c("-10","-6", "-4","-2", "0", "2", "4", "6", "10"), cluster_cols = F)
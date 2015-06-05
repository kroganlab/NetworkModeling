source('R/DB/DBFunctions.R')

CELL_LINE = 'MOUSE_LUNG'
STRAIN = 'H1N1'
LFC = 2
FDR = 0.001

###############
## GET DATASETS

## mRNA ##

CONDITION='.*'
genomics_mouse_H1N1 = getGenomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'genomics_rnaseq_mouse', condition_1 = '.*', condition_2 = '.*')
prot_ph_mouse_H1N1 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ph_sites_mouse', mapping = 'mygene_mouse_1to1', condition_1 = '.*', condition_2 = '.*')
prot_ub_mouse_H1N1 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ub_sites_mouse', mapping = 'mygene_mouse_1to1', condition_1 = '.*', condition_2 = '.*')
metabol_mouse_H1N1 = getMouseMetabolomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, condition_1 = '.*', condition_2 = '.*')
metabol_mouse_H1N1 = getMouseMetabolomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, condition_1 = '.*', condition_2 = '.*', only_known_metabolites = T, mouse_human_orthology_confidence = 1)
mouse_H1N1 = rbind(genomics_mouse_H1N1, prot_ph_mouse_H1N1, prot_ub_mouse_H1N1, metabol_mouse_H1N1)

# library(reshape2)
# mouse_H1N1_w = dcast(entrez_id ~ omics_type+condition_2, data=mouse_H1N1, value.var='log2fc', fun.aggregate = mean)
# mouse_H1N1_w = mouse_H1N1_w[!is.na(mouse_H1N1_w$entrez_id),]
# rownames(mouse_H1N1_w) = mouse_H1N1_w$entrez_id
# mouse_H1N1_w = mouse_H1N1_w[,-1]
# pairs(mouse_H1N1_w)
# library(pheatmap)
# mouse_H1N1_w[is.na(mouse_H1N1_w)]=0
# pheatmap(mouse_H1N1_w)
# pairs(mouse_H1N1_w, upper.panel = panel.cor, pch=21)
# mouse_H1N1_w_c = cor(mouse_H1N1_w, use = 'pairwise.complete.obs', method='pearson')
# pheatmap(mouse_H1N1_w_c)

## get all these from the DB, regardless of cutoff
mouse_H1N1_entrez_unique = unique(mouse_H1N1$entrez_id)
na.omit(mouse_H1N1_entrez_unique)

genomics_mouse_H1N1_consensus = getGenomics(cell_line = CELL_LINE, strain = STRAIN, table = 'genomics_rnaseq_mouse', condition_1 = '.*', condition_2 = '.*', by = "entrez_id", id_list = mouse_H1N1_entrez_unique)
prot_ph_mouse_H1N1_consensus = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, table = 'proteomics_ph_sites_mouse', mapping = 'mygene_mouse_1to1', condition_1 = '.*', condition_2 = '.*', by="entrez_id", id_list = mouse_H1N1_entrez_unique)

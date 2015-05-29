library(data.table)
library(reshape2)
source('R/DB/DBFunctions.R')
source('R/Networks/NetworkFunctions.R')

## CONFIG 
NETWORK_DIR = '~/Projects/HPCKrogan/Data/FluOMICS/database/networks/'
CELL_LINE = 'MOUSE_LUNG'
STRAIN = 'H1N1'
LFC = 2
FDR = 0.05

## GET BASEMAP
mouse_basemap = getMouseBaseMap(basemap_name = 'basemap_menche_mouse', type = '.*', geneid_1_homology_confidence = 1, geneid_2_homology_confidence = 1, geneid_1_human_identity = 0, geneid_2_human_identity = 0)

#################################
## GET DATASETS AND MAKE NETWORKS

## UB ##

CONDITION='H1N1_D00H12'
proteomics_ub_0_12 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ub_sites_mouse', mapping = 'mygene_mouse_1to1', unique_proteins = T, sample_2 = '.*', sample_1 = CONDITION)
writeCytoscapeFile(basemap = mouse_basemap, dataset = proteomics_ub_0_12, filename =sprintf('%smouse_ub_%s.txt',NETWORK_DIR,CONDITION))

CONDITION='H1N1_D01'
proteomics_ub_1_00 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ub_sites_mouse', mapping = 'mygene_mouse_1to1', unique_proteins = T, sample_2 = '.*', sample_1 = CONDITION)
writeCytoscapeFile(basemap = mouse_basemap, dataset = proteomics_ub_1_00, filename =sprintf('%smouse_ub_%s.txt',NETWORK_DIR,CONDITION))

CONDITION='H1N1_D02'
proteomics_ub_2_00 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ub_sites_mouse', mapping = 'mygene_mouse_1to1', unique_proteins = T, sample_2 = '.*', sample_1 = CONDITION)
writeCytoscapeFile(basemap = mouse_basemap, dataset = proteomics_ub_2_00, filename =sprintf('%smouse_ub_%s.txt',NETWORK_DIR,CONDITION))

CONDITION='H1N1_D03'
proteomics_ub_3_00 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ub_sites_mouse', mapping = 'mygene_mouse_1to1', unique_proteins = T, sample_2 = '.*', sample_1 = CONDITION)
writeCytoscapeFile(basemap = mouse_basemap, dataset = proteomics_ub_3_00, filename =sprintf('%smouse_ub_%s.txt',NETWORK_DIR,CONDITION))

CONDITION='H1N1_D04'
proteomics_ub_4_00 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ub_sites_mouse', mapping = 'mygene_mouse_1to1', unique_proteins = T, sample_2 = '.*', sample_1 = CONDITION)
writeCytoscapeFile(basemap = mouse_basemap, dataset = proteomics_ub_4_00, filename =sprintf('%smouse_ub_%s.txt',NETWORK_DIR,CONDITION))

## PH ##

CONDITION='H1N1_D00H12'
proteomics_ph_0_12 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ph_sites_mouse', mapping = 'mygene_mouse_1to1', unique_proteins = T, sample_2 = '.*', sample_1 = CONDITION)
writeCytoscapeFile(basemap = mouse_basemap, dataset = proteomics_ph_0_12, filename =sprintf('%smouse_ph_%s.txt',NETWORK_DIR,CONDITION))

CONDITION='H1N1_D01'
proteomics_ph_1_00 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ph_sites_mouse', mapping = 'mygene_mouse_1to1', unique_proteins = T, sample_2 = '.*', sample_1 = CONDITION)
writeCytoscapeFile(basemap = mouse_basemap, dataset = proteomics_ph_1_00, filename =sprintf('%smouse_ph_%s.txt',NETWORK_DIR,CONDITION))

CONDITION='H1N1_D02'
proteomics_ph_2_00 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ph_sites_mouse', mapping = 'mygene_mouse_1to1', unique_proteins = T, sample_2 = '.*', sample_1 = CONDITION)
writeCytoscapeFile(basemap = mouse_basemap, dataset = proteomics_ph_2_00, filename =sprintf('%smouse_ph_%s.txt',NETWORK_DIR,CONDITION))

CONDITION='H1N1_D03'
proteomics_ph_3_00 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ph_sites_mouse', mapping = 'mygene_mouse_1to1', unique_proteins = T, sample_2 = '.*', sample_1 = CONDITION)
writeCytoscapeFile(basemap = mouse_basemap, dataset = proteomics_ph_3_00, filename =sprintf('%smouse_ph_%s.txt',NETWORK_DIR,CONDITION))

CONDITION='H1N1_D04'
proteomics_ph_4_00 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ph_sites_mouse', mapping = 'mygene_mouse_1to1', unique_proteins = T, sample_2 = '.*', sample_1 = CONDITION)
writeCytoscapeFile(basemap = mouse_basemap, dataset = proteomics_ph_4_00, filename =sprintf('%smouse_ph_%s.txt',NETWORK_DIR,CONDITION))

## mRNA ##

STRAIN = 'H5N1'
CONDITION='H5N1_04D'
genomics_1_4 = getGenomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'genomics_rnaseq_mouse', mapping = 'mygene_mouse_1to1', sample_2 = '.*', sample_1 = CONDITION)
writeCytoscapeFile(basemap = mouse_basemap, dataset = genomics_1_4, filename =sprintf('%smouse_ge_%s.txt',NETWORK_DIR,CONDITION))
CONDITION='H5N1_96H'
proteomics_ph_4_00 = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ph_sites_mouse', mapping = 'mygene_mouse_1to1', sample_2 = '.*', sample_1 = CONDITION)
writeCytoscapeFile(basemap = mouse_basemap, dataset = proteomics_ph_4_00, filename =sprintf('%smouse_ge_%s.txt',NETWORK_DIR,CONDITION))



# ## SPIA test
# allUb = getModProteomics(cell_line = 'MOUSE_LUNG', strain = 'H1N1', lfc = 0, q = 1, table = 'proteomics_ub_sites', mapping = 'mygene_mouse_1to1')
# unique_genes = unique(allUb$entrez_id)
# test = allUb[!is.na(allUb$entrez_id) & allUb$cell_line=='MOUSE_LUNG'& allUb$strain=='H1N1',]
# test = aggregate(log2fc ~ entrez_id, data=test, FUN = max)
# test_vec = test$log2fc
# names(test_vec) = test$entrez_id
# 
# setwd(paste(getwd(),'/scripts/test/',sep=''))
# res=spia(de=test_vec,all=unique_genes,organism="mmu",nB=2000,plots=F,beta=NULL,combine="fisher",verbose=FALSE)
# 
# 

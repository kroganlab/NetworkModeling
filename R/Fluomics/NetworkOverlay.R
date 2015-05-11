library(data.table)
library(reshape2)
source('R/DB/DBFunctions.R')

makeNetworkByUniprot = function(basemap, dataset){
  ds = data.table(dataset)
  bs = data.table(basemap)
  setnames(ds,'protein','uniprot_ac_1')
  step_1 = merge(bs, ds, by='uniprot_ac_1')  
  step_1 = step_1[,sort(colnames(step_1)), with=F]
  setnames(ds,'uniprot_ac_1','uniprot_ac_2')
  step_2 = merge(bs, ds, by='uniprot_ac_2')
  step_2 = step_2[,sort(colnames(step_2)), with=F]
  setnames(ds,'uniprot_ac_2','protein')
  res = rbind(step_1, step_2)
  return(res)
}

makeNetworkByEntrez = function(basemap, dataset){
  ds = data.table(dataset)
  bs = data.table(basemap)
  setnames(ds,'entrez_id','geneid_1')
  step_1 = merge(bs, ds, by='geneid_1')  
  step_1 = step_1[,sort(colnames(step_1)), with=F]
  setnames(ds,'geneid_1','geneid_2')
  step_2 = merge(bs, ds, by='geneid_2')
  step_2 = step_2[,sort(colnames(step_2)), with=F]
  setnames(ds,'geneid_2','entrez_id')
  res = rbind(step_1, step_2)
  return(res)
}

makeAttributesFromNetwork = function(network, dataset){
  nodes = unique(c(network$geneid_1, network$geneid_2))
  node_names = getEntrezNames(values = nodes)
  node_names_add = merge(node_names, dataset[,c('entrez_id','omics_type','condition_1','condition_2','cell_line','strain','log2fc','adj_pvalue')], all.x=T, by='entrez_id')
  ## replace all NA's by numeric values for 
  node_names_add[is.na(node_names_add$log2fc),]$log2fc=0
  node_names_add[is.na(node_names_add$adj_pvalue),]$adj_pvalue=1
  return(node_names_add)
}

## something is wrong with using entrez ids as keys...
writeCytoscapeFile = function(basemap, dataset, filename){
  network = makeNetworkByEntrez(basemap = basemap, dataset = dataset)
  write.table(network, gsub('.txt','-network.txt',filename), eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
  attributes = makeAttributesFromNetwork(network = network, dataset = dataset)
  write.table(attributes, gsub('.txt','-attributes.txt',filename), eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}

writeCytoscapeFile = function(basemap, dataset, filename){
  network = makeNetworkByEntrez(basemap = basemap, dataset = dataset)
  attributes = makeAttributesFromNetwork(network = network, dataset = dataset)
  network$geneid_1 = sprintf('entrez:%s',network$geneid_1)
  network$geneid_2 = sprintf('entrez:%s',network$geneid_2)
  attributes$entrez_id = sprintf('entrez:%s',attributes$entrez_id)
  write.table(network, gsub('.txt','-network.txt',filename), eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
  write.table(attributes, gsub('.txt','-attributes.txt',filename), eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}

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

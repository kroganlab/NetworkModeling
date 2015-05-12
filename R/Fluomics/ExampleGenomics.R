library(data.table)
library(reshape2)
source('R/DB/DBFunctions.R')
source('R/Networks/NetworkFunctions.R')

## don't forget to execute the follwoing in the terminal (opens port forwarding from localhost to higgs server for connection to the DB)
## ssh -L 3307:localhost:3306 [user]@higgs.ucsf.edu

## CONFIG 
NETWORK_DIR = '~/Box Sync/Projects/FluomicsModeling/results/networks/v1/'
CELL_LINE = 'MOUSE_LUNG'
STRAIN = 'H1N1'
LFC = 2
FDR = 0.05

###############
## GET DATASETS

## mRNA ##

CONDITION='H1N1_12H'
genomics_mouse_H1N1_H12 = getGenomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'genomics_rnaseq_mouse', mapping = 'mygene_mouse_1to1', condition_1 = '.*', condition_2 = CONDITION)


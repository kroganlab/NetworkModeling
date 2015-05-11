library(RMySQL)
library(data.table)
library(stringr)

##############
## IMPORT DATA

## DB CONNECTION
## FIRST OPEN TERMINAL AND TYPE:
## ssh -L 3307:localhost:3306 everschueren@higgs.ucsf.edu
con = dbConnect(RMySQL::MySQL(), dbname = 'FluOMICS', host = '127.0.0.1', port = 3307, default.file = '/etc/my.cnf')

## FILE LIST
data_files = read.delim(file = 'scripts/datasets.txt', stringsAsFactors=F)
#importDatasets()

## DISCONNECT
dbDisconnect(con)

##############
## IMPORT MAPPING

## DB CONNECTION
## FIRST OPEN TERMINAL AND TYPE:
## ssh -L 3307:localhost:3306 everschueren@higgs.ucsf.edu
con = dbConnect(RMySQL::MySQL(), dbname = 'FluOMICS', host = '127.0.0.1', port = 3307, default.file = '/etc/my.cnf')

# file = 'database/resources/human_uniprot_mygene.csv'
# res = importProteinToGeneMapping(con=con, file=file, db_table='mygene_human_1to1')

file = 'database/resources/mouse_uniprot_mygene.csv'
res = importProteinToGeneMapping(con=con, file=file, db_table='mygene_mouse_1to1')

dbDisconnect(con)

# hmdb_names = fread('database/resources//hmdb//DFout4.txt')
# hmdb_kegg = fread('database/resources//hmdb//DFout3.txt')
# hmdb_desc = fread('database/resources//hmdb//DFout2_edit.txt')
# hmdb_uniprot = fread('database/resources//hmdb//DFout_edit.txt')
# 
# # hmdb = merge(hmdb_names, hmdb_desc, by = 'hmdb_accesion')
# hmdb = merge(hmdb_names, hmdb_kegg, by = 'hmdb_accesion', all=T)
# length(unique(hmdb$hmdb_accesion))
# write.table(hmdb, file='database/resources//hmdb/hmdb_description.txt', eol='\n', sep='\t', quote=F, row.names=F, col.names=T)





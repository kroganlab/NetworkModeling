#library('random'); 
require(stringi)
setwd("~/Box Sync/Doug/projects/FLU_Networking_code/Cyto_shiny_dev_generic_data_mult_2")

#####
#  Mouse_map_orig <- read.csv(file = "Linked_subdirectories/data/basemap_menche_mouse.csv")
#  source_names <- data.frame(source = unique(Mouse_map_orig$source))
#  save(source_names, file = 'source_names')
load('source_names')



#generic_gene_names_file <- read.delim( file = "Linked_subdirectories/data/humannet_uniprot_keys.txt" )
#generic_gene_names <- unique(generic_gene_names_file$gene_1)

# genomics_mouse nrow = 29655, 
# length(unique(genomics_mouse$entrez_id)) = 4500

entrez_ids = sample(10000:90000, 4500, replace=F)
# symbols <- as.vector(randomStrings(n=4500, len=5, digits=TRUE, upperalpha=TRUE,
#                                    loweralpha=TRUE, unique=TRUE, check=TRUE))

symbols <- stri_rand_strings(n=4500, length=c(4,5,6,7), pattern="[A-Za-z0-9]")


DF_ID_Name <- data.frame(entrez_id = entrez_ids, symbol = symbols)

data1 <- data.frame(entrez_id = sample(as.vector(DF_ID_Name$entrez_id), size = 30000, replace = T ))
data2 <- merge(data1, DF_ID_Name)

# So now we need to create some edges
  
edges1 <- sample(entrez_ids, size = 1000, replace = T)
edges2 <- sample(entrez_ids, size = 1000, replace = T)
 
# this line moved to top for save
# source_names <- data.frame(source = unique(Mouse_map_orig$source))
source_names$int <- 1:50
# sources1 <- sample(as.vector(source_names), size = 1000, replace = T)
sources_int <- data.frame(int = sample(as.vector(1:50), size = 1000, replace = T))
sources_int2 <- data.frame(int2 = sample(as.vector(1:50), size = 1000, replace = T))

sources <- merge(cbind(sources_int, sources_int2), source_names, sort = F)
Species_Map <- data.frame(geneid_1 = as.vector(edges1), geneid_2 = as.vector(edges2), 
                          source = as.vector(sources$source) )

DF <- data.frame(time = sample(c('03H', '01H', '04H', '02H','05H'), size = 30000, replace = T),
                 strain = sample(c('H1N1', 'H3N2', 'H5N1'), size = 30000, replace = T),
                 cell_line = sample(c('HTBE', 'MDM'), size = 30000, replace = T),
                 log2fc = runif(30000, 1.0 , 10.0 ),
                 q_value = runif(30000, .0001 , .1 ))

DF$condition_2 <- paste0(DF$strain, '_' , DF$time)

DF<- cbind(DF, data2)
DF <- DF[ , !names(DF) %in% c("dummy", "time")]
DF$uniprot_ac <- DF$symbol
## set up for write
genomics_human <- DF
Human_map <- Species_Map
#  save(genomics_mouse, file = 'Linked_subdirectories/data/mouse_data_large3')
#  write.table(Mouse_map, file = 'Linked_subdirectories/data/basemap_menche_mouse.csv', row.names = F, sep = "," )
 
 # for human
 
  save(genomics_human, file = 'Linked_subdirectories/data/human_rnaseq_data')
  write.table(Human_map, file = 'Linked_subdirectories/data/basemap_menche.csv', row.names = F, sep = "," )


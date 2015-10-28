
# Read in a file with clusters per gene and determine if there is an overrepresentation in a kegg pathway within a cluster.


#~~~~~~~~~~~~~~~~~~~~
#     LIBRARIES
#~~~~~~~~~~~~~~~~~~~~
suppressMessages(library(reshape2))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(KEGG.db))
source('~/projects/NetworkModeling/R/DB/DBFunctions.R')


#~~~~~~~~~~~~~~~~~~~~
#     FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~
# get number of genes listed in GO for reference
get_ngense_human <- function(){
  entrez_uniprot_map = org.Hs.egUNIPROT
  entrez_uniprot_map_keys = mappedkeys(entrez_uniprot_map)
  entrez_uniprot_keys = as.data.frame(entrez_uniprot_map[entrez_uniprot_map_keys])
  entrez_uniprot_keys_unique = unique(entrez_uniprot_keys$gene_id)  
  ngenes = length(entrez_uniprot_keys_unique)
  return(ngenes)
}


# get the gene symbol and gene name for the matching entrez id's
get_genenames <- function(x){
  #genes <- org.Hs.egGENENAME
  #genes = as.list( genes)
  genes = select(org.Hs.eg.db, keys=as.character(unique(x$entrez_id)),columns=c('SYMBOL', 'GENENAME'), keytype='ENTREZID')
  x = merge(x, genes, by.x='entrez_id', by.y='ENTREZID')
  return(x)
}


# get a list of KEGG id's, terms, and entrez_id's
getKEGGTerms <- function(){
  cat(">> Getting KEGG ID's...\n")
  KEGG = org.Hs.egPATH
  mapped_keys = mappedkeys(KEGG)
  KEGG_Entrez = as.data.frame(KEGG[mapped_keys])
  KEGG = KEGGPATHID2NAME
  mapped_keys = mappedkeys(KEGG)
  KEGG_Terms = as.data.frame(KEGG[mapped_keys])
  KEGG = unique(merge(KEGG_Entrez, KEGG_Terms, by='path_id'))
  
  # DONT NEED THE UNIPROTS
#   # map uniprot_ac into entrez and then to KEGG id's    !!!!!!!!!!!!!!!!!! THIS CAN BE OPTIMIZED!
#   x <- org.Hs.egUNIPROT
#   # Get the entrez gene IDs that are mapped to a Uniprot ID
#   mapped_genes <- mappedkeys(x)
#   # Convert to a list
#   xx <- as.list(x[mapped_genes])
#   cat('>>   Converting uniprot_ac to entrez_id...\n')
#   tmp = melt(xx)
#   names(tmp) = c('uniprot_ac','gene_id')
#   tmp$uniprot_ac = as.character(tmp$uniprot_ac)
#   KEGG = merge(KEGG, tmp, by='gene_id', all.x=T)
  
  #res = unique(merge(go_ids, KEGG, by.x='ENTREZID', by.y='gene_id'))
  #names(res)[9] = "TERM"
  return(KEGG)
}


# use the mouse-human db in the flu db to convert the mouse to human id's
convert2entrez <- function(x, org){
  # connect to mysql db
  con = myConnect()
  
  # convert m.uniprots & m.entrez to h.entrez
  if(org == "mouse"){
    cat(">>   Converting Mouse Uniprots to Entrez...\n")
    # map m.uniprot to m.entrez to h.entrez
    prot.idx = which(grepl('_', x$labels))
    id_list=unique(x$id[prot.idx])
    query = sprintf("select uniprot_ac, entrez_id from mygene_mouse where uniprot_ac in %s", str_c("('",str_c(id_list,collapse = "','"),"')"))  
    res = dbGetQuery(con, query)
    # keep only one entrez per unirpto
    res$entrez_id = gsub(";.*", "", res$entrez_id)
    res = res[-which(res$entrez_id==""),]
    x = merge(x, res, by.x='id', by.y='uniprot_ac', all.x=T)
    prot.idx = which(!is.na(x$entrez_id))
    x$id[prot.idx] = x$entrez_id[prot.idx]
    x$entrez_id = NULL
    
    # Map m.entrez to h.entrez
    id_list=unique(x$id)
    query = sprintf("select entrez_id_mouse, entrez_id_human from mapping_human_mouse_entrez where orthology_confidence = 1 and entrez_id_mouse in %s", str_c("('",str_c(id_list,collapse = "','"),"')"))  
    res = dbGetQuery(con, query)
    x = merge(x, res, by.x='id', by.y='entrez_id_mouse', all.x=T)
    
    # remove anything without a human entrez id
    x$id = x$entrez_id_human
    x$entrez_id_human = NULL
    x = x[-which(is.na(x$id)),]
    
  }else if(org == "human"){  # need to convert h.uniprot to h.entrez
    cat(">>   Converting human Uniprots to Entrez...\n")
    # map h.uniprot to h.entrez
    prot.idx = which(grepl('_', x$labels))
    id_list=unique(x$id[prot.idx])
    query = sprintf("select uniprot_ac, entrez_id from mygene_human where uniprot_ac in %s", str_c("('",str_c(id_list,collapse = "','"),"')"))  
    res = dbGetQuery(con, query)
    # keep only one entrez per unirpto
    res$entrez_id = gsub(";.*", "", res$entrez_id)
    res = res[-which(res$entrez_id==""),]
    x = merge(x, res, by.x='id', by.y='uniprot_ac', all.x=T)
    # copy over found entrez id's
    prot.idx = which(!is.na(x$entrez_id))
    x$id[prot.idx] = x$entrez_id[prot.idx]
    x$entrez_id = NULL
    
  }
  
  # close db connection
  dbDisconnect(con)
  return(x)
}


# Clean up the id's and convert to human id's to match with kegg
clean_ids <- function(x, org){
  cat(">> Cleaning ID's...\n")
  x$dataset = gsub('\\|\\|\\|.*','',x$labels)
  x$id = gsub('.*\\|\\|\\|','',x$labels)
  x$id = gsub('_.*','',x$id)
  
  # convert mouse uniprots to human
  x = convert2entrez(x, org)
  # remove unlabeled metabolites (ALL as of now)
  remove.idx = which(x$dataset == 'MET')
  if(length(remove.idx)>0){
    x = x[-remove.idx,]
  }
  
  return(x)
  }


# checks each cluster to see if it is enriched for a certain term
clusterEnrichment <- function(x, ngenes, kegg){
  cat(">> Checking for term enrichment within each cluster...\n")
  # number of genes in all of kegg
  nkegg = length(unique(kegg$gene_id))
  # number of proteins in each kegg term
  prot.term.kegg = unique(kegg[,c('path_id','gene_id')])
  prot.term.kegg = data.frame(table(prot.term.kegg$path_id))
  
  # itterate through the clusters
  results = c()
  for(i in unique(x$cluster)){
    idx = which(x$cluster == i)
    prot.clust = length(unique(x$entrez_id[idx]))
    prot.term.clust = data.frame(table(x$path_id[idx]))
    prots = merge(prot.term.clust, prot.term.kegg, by='Var1')
    names(prots) = c('path_id','prot.term.clust','prot.term.kegg')
    prots  = cbind(prots, cluster=i, nclust=prot.clust, nkegg=nkegg)
    
    prots$p_val = apply(prots[,-1], 1, function(y){ return( phyper(y[1]-1, y[2], y[5]-y[2], y[4] , lower.tail=F )) })
    results=rbind(results, prots)
  }
  # adjust for multiple comparisons
  results$q_val = p.adjust(results$p_val, method="BH")
  tmp = results[order(results$q_val, decreasing=F),]
  tmp = tmp[,c("cluster", "prot.term.clust","prot.term.kegg", "nclust","nkegg", "p_val", 'q_val', "path_id")]
  return(tmp)
}



cluster_enrichment.main <- function(dat, out_file, org){
  names(dat)[grep('id', names(dat))] = 'labels'
  dat = dat[,c('cluster','labels')]
  
  cat('>> Getting KEGG terms...\n')
  ngenes = get_ngense_human()
  kegg = getKEGGTerms()
  
  # Clean up the id's and convert to human id's to match with kegg
  dat.labeled = clean_ids(dat, org)
  names(dat.labeled)[grep('id', names(dat.labeled))] = 'entrez_id'
  
  # Get kegg id's associated with each entrez_id
  dat.kegg = merge(dat.labeled, kegg, by.x='entrez_id', by.y='gene_id')
  
  # check for pathway enrichment in the clusters
  dat.results = clusterEnrichment(dat.kegg, ngenes, kegg)
  # add keg terms to the id's
  dat.results = merge(dat.results, unique(kegg[,c('path_id','path_name')]), by='path_id')
  
  # get gene names
  cat(">> Annotating pathways...\n")
  dat.kegg = get_genenames(dat.kegg)
  # list gene names in data in each pathway
  gene_names = aggregate(data=dat.kegg, SYMBOL~path_id, FUN=function(y){paste(unique(y), collapse=',')})
  dat.results = merge(dat.results, gene_names, by='path_id')
  
  # get POSSIBLE metabolites in pathway (basically a list of what clusters the diff metabolites show up in)
  dat$dataset = gsub('\\|\\|\\|.*','',dat$labels)
  dat$id = gsub('.*\\|\\|\\|','',dat$labels)
  metabolites = aggregate(data=dat[dat$dataset=='MET',], id~cluster, FUN=function(y){paste(unique(y), collapse=',')})
  names(metabolites)[2] = "possible_metabolite_members"
  dat.results = merge(dat.results, metabolites, by='cluster', all.x=T)
  
  
  # re-order and write out the results
  dat.results = dat.results[order(dat.results$cluster, dat.results$p_val),]
  out_file = gsub('.pdf', '_clusterEnrichment.txt', out_file)
  write.table(dat.results, out_file, quote=F, row.names=F, sep='\t')
  
  cat(">> Enrichment Analysis Complete !!!\n")
  return(dat.results)
}




















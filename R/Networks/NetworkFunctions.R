

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
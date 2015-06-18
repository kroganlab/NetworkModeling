library(RMySQL)
library(data.table)
library(stringr)

myConnect = function(dbname='FluOMICS'){
  con = dbConnect(RMySQL::MySQL(), dbname = dbname, host = '127.0.0.1', port = 3307, default.file = '/etc/my.cnf')
  return(con)
}

getTableNames = function(){
  con = myConnect()
  res = dbGetQuery(con, 'show tables;')
  dbDisconnect(con)
  return(res)
}

getGenomics = function(cell_line, strain, lfc=2, q=0.05, table, condition_1='.*', condition_2='.*', by="confidence", id_list=NULL){
  con = myConnect()
  if(by=="confidence"){
    criteria = sprintf("and R.`q_value` < %s and abs(R.`log2fc`) > %s", q, lfc)
  }else if(by=="entrez_id"){
    criteria = sprintf("and R.entrez_id in %s", str_c("('",str_c(id_list,collapse = "','"),"')"))
  }else{
    criteria = ""
  }
  query = sprintf("select E.experiment_id, E.omics_type, E.condition_1, E.condition_2, E.cell_line, E.strain, E.status, R.entrez_id as 'id', R.entrez_id, R.symbol, '' as uniprot_ac, '' as 'kegg_id', R.log2fc, R.q_value as 'q_value' from experiments E join `%s` R on R.`experiment_id` = E.`experiment_id` and R.`sample_1` = E.`sample_1_code` and R.`sample_2` = E.`sample_2_code` where E.`cell_line` regexp '%s' and E.`strain` regexp '%s' and E.`condition_1` regexp '%s' and E.`condition_2` regexp '%s' %s", table, cell_line, strain, condition_1, condition_2, criteria)
  #cat(str_c(query,'\n'))
  res = dbGetQuery(con, query)
  dbDisconnect(con)
  return(res)
}

getModProteomics = function(cell_line, strain, lfc=2, q=0.05, table, mapping='mygene_human', condition_1='.*', condition_2='.*', by="confidence", id_list=NULL){
  con = myConnect()
  if(by=="confidence"){
    criteria = sprintf("and P.`adj_pvalue` < %s and abs(P.`log2fc`) > %s", q, lfc)
  }else if(by=="entrez_id"){
    criteria = sprintf("and GP.entrez_id in %s", str_c("('",str_c(id_list,collapse = "','"),"')"))
  }else if(by=="id"){
    criteria = sprintf("and P.mod_sites in %s", str_c("('",str_c(id_list,collapse = "','"),"')"))
  }else{
    criteria = ""
  }
  query = sprintf("select E.experiment_id, E.omics_type, E.condition_1, E.condition_2, E.cell_line, E.strain, E.status, P.mod_sites as 'id', GP.entrez_id, P.symbol, P.protein as uniprot_ac, '' as 'kegg_id', P.log2fc, P.adj_pvalue as 'q_value' from experiments E join `%s` P on P.`experiment_id` = E.`experiment_id` and P.`sample_1` = E.`sample_1_code` and P.`sample_2` = E.`sample_2_code` left join %s GP on GP.`uniprot_ac` = P.`protein` where E.`cell_line` regexp '%s' and E.`strain` regexp '%s' and E.`condition_1` regexp '%s' and E.`condition_2` regexp '%s' %s", table, mapping, cell_line, strain, condition_1, condition_2, criteria)
  #cat(str_c(query,'\n'))
  res = dbGetQuery(con, query)
  dbDisconnect(con)
  return(res)
}

# special function for mouse because mapping is complicated 
getMouseMetabolomics = function(cell_line, strain, lfc, q, condition_1='.*', condition_2='.*', only_known_metabolites=FALSE, mouse_human_orthology_confidence=1, by="confidence", id_list=NULL){
  con = myConnect()
  if(by=="confidence"){
    criteria = sprintf("and M.`adj_pvalue` < %s and abs(M.`log2fc`) > %s", q, lfc)
  }else if(by=="entrez_id"){
    criteria = sprintf("and HE.`entrez_gene_id_mouse` in %s", str_c("('",str_c(id_list,collapse = "','"),"')"))
  }else if(by=="id"){
    criteria = sprintf("and M.metabolite in %s", str_c("('",str_c(id_list,collapse = "','"),"')"))
  }else{
    criteria = ""
  }
  if(only_known_metabolites){
    query = sprintf("select E.experiment_id, E.omics_type, E.condition_1, E.condition_2, E.cell_line, E.strain, E.status, M.metabolite as 'id', M.rounded_mass_id, HE.`entrez_gene_id_mouse` as 'entrez_id', '' as 'uniprot_ac', MI.`kegg_id`, M.log2fc, M.adj_pvalue as 'q_value' from experiments E join `metabolomics_mouse` M on M.`experiment_id` = E.`experiment_id` and M.`sample_1` = E.`sample_1_code` and M.`sample_2` = E.`sample_2_code` join `metabolomics_ids` MI on MI.`rounded_mass_id` = M.`rounded_mass_id` join `mapping_hmdb_mouse_entrez` HE on HE.`hmdb` = MI.`hmdb_id`  where E.`cell_line` regexp '%s' and E.`strain` regexp '%s' and E.`condition_1` regexp '%s' and E.`condition_2` regexp '%s' and HE.`orthology_confidence` >= %s %s", cell_line, strain, condition_1, condition_2, mouse_human_orthology_confidence, criteria)   
  }else if(only_known_metabolites==F){
    query = sprintf("select E.experiment_id, E.omics_type, E.condition_1, E.condition_2, E.cell_line, E.strain, E.status, M.metabolite as 'id', M.rounded_mass_id, HE.`entrez_gene_id_mouse` as 'entrez_id', '' as 'uniprot_ac', MI.`kegg_id`, M.log2fc, M.adj_pvalue as 'q_value' from experiments E join `metabolomics_mouse` M on M.`experiment_id` = E.`experiment_id` and M.`sample_1` = E.`sample_1_code` and M.`sample_2` = E.`sample_2_code` left join `metabolomics_ids` MI on MI.`rounded_mass_id` = M.`rounded_mass_id` left join `mapping_hmdb_mouse_entrez` HE on HE.`hmdb` = MI.`hmdb_id`  where E.`cell_line` regexp '%s' and E.`strain` regexp '%s' and E.`condition_1` regexp '%s' and E.`condition_2` regexp '%s' %s", cell_line, strain, condition_1, condition_2, criteria)    
  }else{
    cat('bad combination of parameters\n')
    quit()
  }
  cat(str_c(query,'\n'))
  res = dbGetQuery(con, query)
  dbDisconnect(con)
  return(res)
}

## WARNING : SO FAR WE'RE NOT MAPPING TO ANYTHING YET, WE'RE JUST USING THE PARAM AS A PLACEHOLDER FOR LATER
# getMetabolomics = function(cell_line, strain, lfc, q, table, mapping='mapping_hmdb_human_uniprot', condition_1='.*', condition_2='.*'){
#   con = myConnect()
#   query = sprintf("select E.experiment_id, E.omics_type, E.condition_1, E.condition_2, E.cell_line, E.strain, E.status,'' as 'kegg_id', M.metabolite, '' as entrez_id, M.log2fc, M.adj_pvalue as 'q_value' from experiments E join `%s` M on M.`experiment_id` = E.`experiment_id` and M.`sample_1` = E.`sample_1_code` and M.`sample_2` = E.`sample_2_code` where E.`cell_line` regexp '%s' and E.`strain` regexp '%s' and M.`adj_pvalue` < %s and abs(M.`log2fc`) > %s and E.`condition_1` regexp '%s' and E.`condition_2` regexp '%s'", table, cell_line, strain, q, lfc, condition_1, condition_2)
#   cat(str_c(query,'\n'))
#   res = dbGetQuery(con, query)
#   dbDisconnect(con)
#   return(res)
# }


getBaseMap = function(basemap_name){
  con = myConnect()
  query = sprintf("select B.*, MH1.`uniprot_ac` as 'uniprot_ac_1', MH2.`uniprot_ac` as 'uniprot_ac_2' from %s B join `mygene_human_1to1` MH1  on MH1.`entrez_id` = B.`geneid_1` join `mygene_human_1to1` MH2 on MH2.`entrez_id` = B.`geneid_2` where MH1.`uniprot_ac_reviewed` = TRUE and MH2.`uniprot_ac_reviewed` = TRUE",basemap_name)
  res = dbGetQuery(con, query)
  dbDisconnect(con)
  return(res)
}

getNames = function(name_table, values, original_id, new_id='gene_name'){
  con = myConnect()
  query = sprintf("select %s, %s from %s where %s in ('%s') and status='reviewed'",original_id, new_id, name_table, original_id, paste(values, collapse="','"))
  res = dbGetQuery(con, query)
  dbDisconnect(con)
  return(res)
}

writeRNAseqGenomics = function(con, data, db_table){
  data_str = paste(apply(data,1,function(x)paste(x,collapse="','")),collapse="'),('")
  query = paste(sprintf("INSERT INTO %s VALUES ('%s')", db_table, data_str))
  res = dbGetQuery(con, query)
  cat(sprintf('INSERTED %s rows for EXPERIMENT %s in TABLE %s \n', nrow(data), unique(data$experiment_id), db_table))
}

writeProject3Genomics = function(con, data, db_table){
  data[,test_id:=NULL]
  data_str = paste(apply(data,1,function(x)paste(x,collapse="','")),collapse="'),('")
  query = paste(sprintf("INSERT INTO %s VALUES ('%s')", db_table, data_str))
  res = dbGetQuery(con, query)
  cat(sprintf('INSERTED %s rows for EXPERIMENT %s in TABLE %s \n', nrow(data), unique(data$experiment_id), db_table))
}

writePhosphoProteomics = function(con, data, db_table){
  data=data[,c('Protein','mod_sites','log2FC','SE','Tvalue','DF','pvalue','adj.pvalue','entrezgene','uniprot_genename','sample_1','sample_2','experiment_id'),with=F]
  data_str = paste(apply(data,1,function(x)paste(x,collapse="','")),collapse="'),('")
  query = paste(sprintf("INSERT INTO %s VALUES ('%s')", db_table, data_str))
  res = dbGetQuery(con, query)
  cat(sprintf('INSERTED %s rows for EXPERIMENT %s in TABLE %s \n', nrow(data), unique(data$experiment_id), db_table))
}

writeUbProteomics = function(con, data, db_table){
  data=data[,c('Protein','mod_sites','log2FC','SE','Tvalue','DF','pvalue','adj.pvalue','entrezgene','uniprot_genename','sample_1','sample_2','experiment_id'),with=F]
  data_str = paste(apply(data,1,function(x)paste(x,collapse="','")),collapse="'),('")
  query = paste(sprintf("INSERT INTO %s VALUES ('%s')", db_table, data_str))
  res = dbGetQuery(con, query)
  cat(sprintf('INSERTED %s rows for EXPERIMENT %s in TABLE %s \n', nrow(data), unique(data$experiment_id), db_table))
}

writeAPMSProteomics = function(con, data, db_table){
  data=data[,c('Protein','log2FC','SE','Tvalue','DF','pvalue','adj.pvalue','entrezgene','uniprot_genename','sample_1','sample_2','experiment_id'),with=F]
  data_str = paste(apply(data,1,function(x)paste(x,collapse="','")),collapse="'),('")
  query = paste(sprintf("INSERT INTO %s VALUES ('%s')", db_table, data_str))
  res = dbGetQuery(con, query)
  cat(sprintf('INSERTED %s rows for EXPERIMENT %s in TABLE %s \n', nrow(data), unique(data$experiment_id), db_table))
}

writeMetabalomics = function(con, data, db_table){
  experiment_id = data$experiment_id
  data=data[,2:19,with=F]
  data[,sample_1:='']
  data[,sample_2:='']
  data[,experiment_id:=experiment_id]
  data_str = paste(apply(data,1,function(x)paste(dbEscapeStrings(con,x),collapse="','")),collapse="'),('")
  query = paste(sprintf("INSERT INTO %s VALUES ('%s')", db_table, data_str))
  res = dbGetQuery(con, query)
  cat(sprintf('INSERTED %s rows for EXPERIMENT %s in TABLE %s \n', nrow(data), unique(data$experiment_id), db_table))
}

importProteinToGeneMapping = function(con, file, db_table){
  data = fread(file, sep = '\t')
  entrezids = c()
  uniprotacs = c()
  for(i in 1:nrow(data)){
    entry = data[i,]
    splitlist = as.numeric(unlist(str_split(entry$entrezgene, ';')))
    entrezids = c(entrezids, splitlist)
    uniprotacs = c(uniprotacs, rep(entry$Entry,length(splitlist)))
  }
  res = data.table(uniprotacs=uniprotacs,entrezids=entrezids)
  res = na.omit(res)
  data_str = paste(apply(res,1,function(x)paste(x,collapse="','")),collapse="'),('")
  query = paste(sprintf("INSERT INTO %s VALUES ('%s')", db_table, data_str))
  dbGetQuery(con, query)
  cat(sprintf('INSERTED %s rows for MAPPING %s \n', nrow(res), db_table))
}

importDatasets = function(con){
  ## genomics data
  for(i in 1:nrow(data_files)){
    entry = data_files[i,]
    data = fread(entry$file)
    data[,significant:=NULL]
   
    filename = gsub('\\.(csv|txt|diff)','',last(unlist(str_split(entry$file, '/'))),)
    data[,experiment_id:=filename]
   
    if(entry$import == 1){
      cat(sprintf("READ %s\n",entry$file))
      
      if(entry$table == 'genomics_rnaseq'){
        writeRNAseqGenomics(con, data, entry$table)
      }else if(entry$table == 'proteomics_ph_sites'){
        writePhosphoProteomics(con, data, entry$table)
      }else if(entry$table == 'proteomics_ub_sites'){
        writeUbProteomics(con, data, entry$table)
      }else if(entry$table == 'proteomics_apms'){
        writeAPMSProteomics(con, data, entry$table)
      }else if(entry$table == 'genomics_rnaseq_p3'){
        writeProject3Genomics(con, data, entry$table)
      }else if(entry$table == 'metabolomics'){
        writeMetabalomics(con, data, entry$table)
      }
      
      tryCatch({
        query = sprintf("INSERT INTO data_import VALUES ('%s','%s','%s')", dbEscapeStrings(con, entry$file), entry$table, 1)
        res = dbGetQuery(con, query)
      },
      error=function(cond) {
        message(paste("ENTRY already exists:\n", entry$file))
        message(cond)
        # Choose a return value in case of error
        return(NULL)
      },
      warning=function(cond) {
        message(paste("ENTRY already exists:\n", entry$file))
        message(cond)
        # Choose a return value in case of error
        return(NULL)
      },
      finally={
        
      })
      
    }
  }
}


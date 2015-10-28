# This script requires the user to have forwarded their MYSQL port 
source('~/projects/NetworkModeling/R/DB/DBFunctions.R')
source('~/projects/NetworkModeling/R/Networks/NetworkFunctions.R')
source('~/projects/NetworkModeling/R/Fluomics/clustering/top_down_dend_clustering.R')
source('~/projects/NetworkModeling/R/Fluomics/clustering/clustergrams.R')

library(ggplot2)
library(data.table)
library(reshape2)
library(RMySQL)
library(stringr)
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))


# ~~~~~~~~~~ BEGIN FUNCTIONS ~~~~~~~~~~

# Grab ALL time points involved with a gene with ANY significant time point
globalCluster.getDatasets <- function(CELL_LINE, STRAIN, LFC=2, FDR=0.001, CONDITION='.*'){
  
  if(grepl('MOUSE', CELL_LINE)){
    # get all the significant entrez id's from the sets
    genomics_mouse = getGenomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'genomics_rnaseq_mouse', condition_1 = '.*', condition_2 = '.*', by = "confidence")
    prot_ph_mouse = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ph_sites_mouse', mapping = 'mygene_mouse_1to1', condition_1 = '.*', condition_2 = '.*', by = "confidence")
    prot_ub_mouse = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, table = 'proteomics_ub_sites_mouse', mapping = 'mygene_mouse_1to1', condition_1 = '.*', condition_2 = '.*', by = "confidence")
    metabol_mouse_known_only = getMouseMetabolomics(cell_line = CELL_LINE, strain = STRAIN, lfc = LFC, q = FDR, condition_1 = '.*', condition_2 = '.*', only_known_metabolites = F, mouse_human_orthology_confidence = 1)
    #mouse_all = rbind(genomics_mouse, prot_ph_mouse, prot_ub_mouse, metabol_mouse_known_only)
    
    
    # use the specific unique identifier for each set instead of entrez in order to avoid pulling multiple sites w entrez
    # pull all time points for any gene with a significant time point
    genomics_mouse_consensus = getGenomics(cell_line = CELL_LINE, strain = STRAIN, table = 'genomics_rnaseq_mouse', condition_1 = '.*', condition_2 = '.*', by = "entrez_id", id_list = na.omit(unique(genomics_mouse$entrez_id)) )
    prot_ph_mouse_consensus = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, table = 'proteomics_ph_sites_mouse', mapping = 'mygene_mouse_1to1', condition_1 = '.*', condition_2 = '.*', by="id", id_list = na.omit(unique(prot_ph_mouse$id)) )
    prot_ub_mouse_consensus = getModProteomics(cell_line = CELL_LINE, strain = STRAIN, table = 'proteomics_ub_sites_mouse', mapping = 'mygene_mouse_1to1', condition_1 = '.*', condition_2 = '.*', by="id", id_list = na.omit(unique(prot_ub_mouse$id)) )
    # Metabolomics!
    prot_mb_mouse_consensus = getMouseMetabolomics(cell_line = CELL_LINE, strain = STRAIN, condition_1 = '.*', condition_2 = '.*', by="id", id_list = na.omit(unique(metabol_mouse_known_only$id)), only_known_metabolites = F)
    # remove the "duplicates" created when joining with the metabolomics_ids table to get entrez_id's 
    prot_mb_mouse_consensus$pairs = paste(prot_mb_mouse_consensus$id, prot_mb_mouse_consensus$condition_2,sep='|')
    prot_mb_mouse_consensus = prot_mb_mouse_consensus[!duplicated(prot_mb_mouse_consensus$pairs),]
    prot_mb_mouse_consensus$pairs = prot_mb_mouse_consensus$rounded_mass_id = c()
    prot_mb_mouse_consensus$symbol = prot_mb_mouse_consensus$id
    
    all_consensus = rbind(genomics_mouse_consensus, prot_ph_mouse_consensus, prot_ub_mouse_consensus, prot_mb_mouse_consensus)
  }else{
  
    ## UPDATE FOR HUMAN
    
  }
  
  return(all_consensus)
}


# creates a column name with cluster id and the percentages of each cluster made up by each dataset
clusterBreakdown <- function(x){
  cat('Calculating dataset contributions to each cluster...\n')
  x$dataset = gsub('\\|\\|\\|.*','',x$id)
  
  tmp = dcast(data=x, cluster~dataset)
  tmp[,-1] = round(tmp[,-1]/apply(tmp[,-1], 1, sum),2)*100
  tmp$clustertitle = tmp$id
  # add dataset names to clusters row by row
  for(i in 1:dim(tmp)[1]){ 
    dat = tmp[i,]
    dat[grep('^MET', names(dat))] = paste('MET:', dat[grep('^MET', names(dat))],'%', sep="")
    dat[grep('^PH', names(dat))] = paste('PH:', dat[grep('^PH', names(dat))],'%', sep="")
    dat[grep('^UB', names(dat))] = paste('UB:', dat[grep('^UB', names(dat))],'%', sep="")
    dat[grep('RNAseq', names(dat))] = paste('RNAseq:', dat[grep('RNAseq', names(dat))],'%', sep="")
    
    data.percentages = paste( dat[, -c(grep(":0%", dat), grep('cluster', names(dat)))] , collapse=', '   )
    tmp$clustertitle[i] = paste(dat$cluster, " - ", data.percentages, sep="")
  }
  
  tmp = merge(x, tmp[,c('cluster','clustertitle')], by='cluster')
  return(tmp)
}


# Used to plot median function in the cluster plots
stat_sum_single <- function(fun, geom="point", ...) {
  stat_summary(fun.y=fun, colour="red", geom=geom, size = 1, ...)
}


# creates plots like the STEM cluster plots with the addition of a Median line (dotted red) and the population distribution of the datasets in each cluster (the title)
plotClusters <- function(x, out_file){
  x$id = row.names(x)
  
  # add % of each dataset to cluster name
  x = clusterBreakdown(x)
  
  tmp = melt(x, id=c('cluster', 'clustertitle','id','dataset'), variable.name='time', value.name='log2FC')
  tot_cols = 3
  p = ggplot(tmp, aes(x=time, y=log2FC, group=id) ) + geom_line(aes(colour=dataset)) + facet_wrap(facets= ~clustertitle, scales='free_x', ncol=tot_cols) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_sum_single(median, geom="line", aes(group=cluster, size=1), linetype="dashed") 
  cat(paste('Creating cluster plot at ', out_file, '\n', sep=""))
  ggsave(filename=out_file, plot=p, limitsize=FALSE, width=11, height=11/tot_cols*ceiling(length(unique(tmp$cluster))/tot_cols) )
  
  n=dim(x)[2]
  x=x[,c(n, 1:(n-1))]
  x$dataset = c()
  return(x)
}


# adds a column of gene names to the data
annotateEntrez <- function(x, CELL_LINE){
  if(grepl('MOUSE',CELL_LINE)){
    annotation_table = 'mygene_mouse'
  }else{
    annotation_table = 'mygene_human'
  }
  
  # get entrez
  x$entrez_id_mygene = gsub('[a-zA-Z|_]+-','',x$id)
  # keep only the entries that are in our set
  entrez_id_mygene = unique(unlist(strsplit(x$entrez_id_mygene,',')))
  # get the gene name listing from the db
  con = myConnect()
  res = dbGetQuery(conn = con, statement = sprintf("select entrez_id_mygene, gene_name from %s where entrez_id_mygene regexp '%s' group by entrez_id_mygene, gene_name", annotation_table, paste(entrez_id_mygene, collapse='|')))
  dbDisconnect(con)
  
  # Merge singles
  x = merge(x, res, by ='entrez_id_mygene', all.x=TRUE)
  
  # merge multiples
  idx = grep(',', x$entrez_id_mygene)
  for(i in 1:length(idx)){
    ids = unlist(strsplit(x[idx[i],'entrez_id_mygene'], ','))
    #sum(res$entrez_id_mygene %in% ids)
    genes = c()
    for(j in 1:length(ids)){
      tmp = res[grep(ids[j], res$entrez_id_mygene),]
      genes = c(genes, tmp$gene_name[1])
    }
    genes = paste(unique(genes), collapse=',')
    x[idx[i],'gene_name'] = genes
  }
  
  x$entrez_id_mygene = c()
  return(x)  
}


# removes clusters with only one data type and clusters that have no FC variance (straight lines)
cleanClusters <- function(x, out_file, CELL_LINE){
  thresh.diff = 1
  tmp=x
  tmp$id = row.names(x)
  
  cat('Filtering clusters...\n')
  # per cluster: calculate the median values for each time point, then check the variance between them all
  tmp = melt(data=tmp, id=c('id','cluster'), variable.name='time', value.name='log2FC')
  tmp = dcast(tmp[,-1], cluster~time, value.var='log2FC', median)
  tmp = as.data.frame( cbind(tmp, diff=apply(tmp[,-1], 1, function(y){return(max(abs(y)))})), stringsAsFactors=F)
  
  idx = which(tmp$diff > thresh.diff)  # index of where variance meets threshold
  x.new = x[x$cluster %in% tmp[idx,'cluster'],]  # keep only clusters that meet the threshold
  x.remove = x[x$cluster %in% tmp[-idx,'cluster'],]  # keep only clusters that FAIL to meet the threshold
  
  # Plot and save the clusters that meet the threshold
  cat('Saving clusters meeting the threshold...\n')
  x.new = plotClusters( x.new, out_file )
  # Plot and save the uniportant "Removed" clusters
  if(dim(x.remove)[1]>0){
    cat('Saving clusters failing to meet the threshold...\n')
    x.remove = plotClusters( x.remove, gsub('.pdf','_REMOVED.pdf', out_file) )
    
    # combine sets
    x.new$threshold = 1
    x.remove$threshold = 0
    x.all = rbind(x.new, x.remove)
  }else{
    x.all = x.new
  }
  ###cat('Annotating data file...\n')
  ###x.all = annotateEntrez(x.all, CELL_LINE)
  
  write.table(x.all, gsub('.pdf','_data.txt',out_file), quote=F, row.names=F, col.names=T, sep='\t')
  return(x.all)
}


# Remove the metabolite and entries that had a one to many mapping to entrez id's so we only count the single entry in the clustering
fixIDs <- function(x){
  #Shorten omics type
  x$omics_type[grep('rnaseq', x$omics_type)] = 'RNAseq'
  x$omics_type[grep('ph', x$omics_type)] = 'PH'
  x$omics_type[grep('ub', x$omics_type)] = 'UB'
  x$omics_type[grep('metabolomics', x$omics_type)] = 'MET'
  
  x$id = paste(x$omics_type,x$id, sep="|||")
  tmp = aggregate(data=x[,-1], id~., function(y){paste(unique(y), collapse=',')} )
  tmp$id = gsub(',[a-zA-z]+\\|\\|\\|',',', tmp$id)
  return(tmp)
}


globalCluster.main <- function(dat, out_file, cluster_size=50, cluster_method='kmeans', cor.thresh=.9){
  dat.wide = dcast(data=dat, omics_type+id~condition_2, value.var='log2fc', median, na.rm=T)
  
  cat('Removing one to many conversion duplicates...\n')
  # combine all entrez id's of the one to many mappings so we can remove the duplicate values for clustering.
  dat.wide = fixIDs(dat.wide)
  
  # set rownames as identifiers
  row.names(dat.wide) = dat.wide$id
  dat.wide$id = c()
  
  ## CLUSTER DATA
  #--------------------------------------
  if(cluster_method == 'kmeans'){
    # K-MEANS
    #  Find best number of clusters  --> Found to be 10, but they don't look that great
    #library(NbClust)
    #set.seed(456456)
    #nc <- NbClust(dat, min.nc=10, max.nc=80, method='kmeans', distance="minkowski", index="all")
    #barplot(table(nc$Best.n[1,]), xlab="Numer of Clusters", ylab="Number of Criteria", main="Number of Clusters Chosen by 26 Criteria")
    cat(paste('Clustering using ', cluster_size, ' clusters...\n', sep=""))
    set.seed(456456)
    fit = kmeans(dat.wide, centers=cluster_size, iter.max = 100, nstart=10)
    dat.wide = cbind(dat.wide, cluster=fit$cluster)
    clusters.good = cleanClusters(dat.wide,out_file, CELL_LINE)
    
    #clustergram  --  to help determine how many clusters to us
    cat('Creating clustergram...\n')
    tmp <- scale(dat.wide)
    pdf(gsub('.pdf','_clustergram.pdf', out_file))
    clustergram(tmp, k.range = 1:50, line.width = 0.004) # notice how I am using line.width.  Play with it on your problem, according to the scale of Y.
    dev.off()
    
  }else if(cluster_method=='topdown'){
    # perform top-down clustering on the dendocram created from hierarchical clustering based on a correlation matrix.
    clusters.good = topdown.main(dat.wide, cor.thresh=cor.thresh, data_file=out_file)
    
  }
  cat('Clustering Complete!\n')
  return(clusters.good)
  
}



























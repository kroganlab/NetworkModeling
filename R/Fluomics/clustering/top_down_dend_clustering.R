#! /usr/bin/Rscript --vanilla
suppressMessages(library(dendextend))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))

#~~~~~~~~~~~~~~~~~~~~~~~~
#     FUNCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~

# Traverse a dendogram based on a correlation matrix and make clusters based on if the children have a correlation within the threshold
dend.walk <- function(x.dend, x.cor, cor.thresh){
  #print(attributes(x.dend)$members)
  cat('-')
  
  # get location/height to cut the tree at to separate into 2 parts
  dend.height = sort(x.dend %>% hang.dendrogram %>% get_nodes_attr("height"), decreasing=T)[-1]
  
  # cut the tree
  tmp.tree = cutree(x.dend, h=dend.height[1])
  
  # get correlation of cluster 1 members
  idx1 = names(tmp.tree)[tmp.tree == 1]
  if( min(x.cor[idx1, idx1]) < cor.thresh ){    # if all of the children have a low correlation, then cut them again
    x1.dend = prune(x.dend, leaves=c(names(tmp.tree)[tmp.tree != 1]))
    tmp1 = dend.walk(x1.dend, x.cor, cor.thresh)
  }else{  # use unique number as a cluster number
    cat('-\n')
    tmp1 = cbind(cluster = runif(1, 0, 1), idx1)
  }
  
  # get correlation of cluster 2 members
  idx2 = names(tmp.tree)[tmp.tree == 2]
  if( min(x.cor[idx2, idx2]) < cor.thresh ){    # if all of the children have a low correlation, then cut them again
    x2.dend = prune(x.dend, leaves=c(names(tmp.tree)[tmp.tree != 2]))
    tmp2 = dend.walk(x2.dend, x.cor, cor.thresh)
  }else{  # use unique number as a cluster number
    cat('-\n')
    tmp2 = cbind(cluster = runif(1, 0, 1), idx2)
  }
  
  clusters.all = rbind(tmp1, tmp2)
  return(clusters.all)
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


# plot out the different proteins/etc in clusters
traverse.plot <- function(tclusts, dat, outfile){
  dat = cbind(id=row.names(dat), dat)
  dat = melt(dat, by="id", variable.name='time', value.name='log2FC')
  dat = merge(dat, tclusts, by.x='id', by.y='labels')
  dat$dataset = gsub('\\|\\|\\|.*','',dat$id)
  
  dat = clusterBreakdown(dat)
  
  cat('Plotting...\n')
  outfile = gsub('.pdf', '_topdown.pdf', outfile)
  print(outfile)
  cat(paste('  Saving plot to: ', outfile, '\n', sep=''))
  tot_cols = 3
  p = ggplot(dat, aes(x=time, y=log2FC, group=id) ) + geom_line(aes(colour=dataset)) + facet_wrap(facets= ~clustertitle, scales='free_x', ncol=tot_cols) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_sum_single(median, geom="line", aes(group=cluster, size=1), linetype="dashed") 
  #cat(paste('Creating cluster plot at ', out_file, '\n', sep=""))
  ggsave(filename=outfile, plot=p, limitsize=FALSE, width=11, height=11/tot_cols*ceiling(length(unique(dat$cluster))/tot_cols) )
}



topdown.main <- function(x, cor.thresh=.9, data_file){
  cat('>> Beginning Hierarchical top down clustering...\n')
  # remove samples that don't change over time:
  idx = which(apply(x,1,var)==0)
  if(length(idx)>0){
    x = x[-idx,]
  }
  
  cat('>> Creating correlation matrix...\n')
  # use correlations in place of distances
  x.cor = cor(t(x))
  print(data_file)
  outfile=gsub('.pdf',paste('_x.cor_',cor.thresh,'.txt', sep=''), data_file)
  write.table(x.cor, outfile, quote=F, row.names=T, col.names=T, sep='\t')
  x.dist = as.dist(1-x.cor, diag=T, upper=F)
  x.clust = hclust(x.dist)
  
  x.dend = as.dendrogram(x.clust)
  #plot(x.dend)
  
  cat('>> Clustering the dendrogram from the top down...\n')
  # get the heights of the tree, make a cut
  tmp = dend.walk(x.dend, x.cor, cor.thresh)
  tmp = as.data.frame(tmp, stringsAsFactors=F)
  tmp$cluster = as.numeric(factor(tmp$cluster, levels=sort(unique(tmp$cluster))))
  names(tmp) = c('cluster','labels')
  
  # write out results
  cat('>> Writing out results...\n')
  #outfile=gsub('.pdf', paste( '_top_down_hclust_',cor.thresh,'.txt', sep=''), data_file)
  outfile = paste(dirname(data_file), '/top_down_hclust_',cor.thresh,'.txt', sep='')
  write.table(tmp, outfile, quote=F, row.names=T, col.names=T, sep='\t')
  cat('>> Plotting results...\n')
  temp = traverse.plot(tmp, x, gsub('.txt','.pdf', outfile))
  return(tmp)
}


















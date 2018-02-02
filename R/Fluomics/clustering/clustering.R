library(ggplot2)
library(reshape2)
library(plotly)

# setwd('~/Box Sync/projects/cardiomyocytes/results')
setwd('~/Box Sync/projects/FluomicsModeling/data/datasets-fluomics/metabolomics_colaboration/')   # <-------------------------------------- INPUT

# what set do you want to run?
out_file <- 'plots/Clustering/kmeans/20170119_kmeans_MOUSE.pdf'                 # <-------------------------------------- INPUT
data_file <- 'data/fluomics_db_20170119.txt'                                    # <-------------------------------------- INPUT
species = 'MOUSE'  # MOUSE|HTBE  -- whatever the "Model" value is               # <-------------------------------------- INPUT

dat <- read.delim(data_file, stringsAsFactors=F)

# fix the names
names(dat)[grep("Log2FC", names(dat))] <- 'log2FC'
names(dat)[grep("Pvalue", names(dat))] <- 'adj.pvalue'
names(dat)[grep("Experiment", names(dat))] <- 'omics_type'

# set up the timepoints to use
dat$condition_2 <- gsub("(.*_)([A-Z0-9]+)(-.*)","\\2", dat$Comparison)

# add the unique identifiers
dat$id <- paste(dat$omics_type, dat$Model, dat$FLU, dat$Protein, dat$Gene, sep='_')
dat.old <- dat




#thresholds
p_val_thresh <- 0.05          # <-------------------------------------- INPUT
log2FC_thresh <- 0.58         # <-------------------------------------- INPUT

# Apply threshold
dat <- dat[ which( (dat$adj.pvalue < p_val_thresh) ),]
dat <- dat[ which( (dat$log2FC > log2FC_thresh) | (dat$log2FC < -log2FC_thresh)),]

# Get the rest of the timepoint data that correslponds with the proteins
prots <- unique(dat$id)
dat <- dat.old[dat.old$id %in% prots,]

### SELECT SPECIES/DATA TYPE
unique(dat$Model)
dat <- dat[which(dat$Model == species ),]




# NOTICE THIS PREVIOUS VERSION!
# source('~/projects/NetworkModeling-preCrash/R/Fluomics/clustering/globalDataClustering.R')
source('~/github/kroganlab/NetworkModeling/R/Fluomics/clustering/globalDataClustering.R')

CLUSTER_SIZE <- 100           # <-------------------------------------- INPUT
CLUSTER_METHOD <- 'kmeans'    # <-------------------------------------- INPUT
COR_THRESH <- 0.9             # <-------------------------------------- INPUT
dat.clusters = globalCluster.main(dat, out_file=out_file, cluster_size=CLUSTER_SIZE, cluster_method=CLUSTER_METHOD, cor_thresh=COR_THRESH)



# create a plot
cluster_of_interest = 46
x <- dat.clusters[which(dat.clusters$cluster==cluster_of_interest),]
x$dataset = gsub("([A-Za-z]+)(_.*)","\\1", x$id)
tmp = melt(x, id=c('cluster', 'clustertitle','id', 'dataset'), variable.name='time', value.name='log2FC')
p = ggplot(tmp, aes(x=time, y=log2FC, group=id) ) + geom_line(aes(colour=dataset)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + stat_sum_single(median, geom="line", aes(group=cluster, size=1), linetype="dashed") 

ggplotly(p)









#~~~~~~~~~~~~~~~~~~~~~~~
# Enrichment Analysis
#~~~~~~~~~~~~~~~~~~~~~~~
source('~/projects/mist/src/Enrichment_v2.R')

# Clean up uniprots
dat.clusters$prey <- gsub('[A-Z]+\\|\\|\\|','',dat.clusters$id)
dat.clusters$prey <- gsub('_.*$','',dat.clusters$prey)

# Create directory to save results in
out_dir <- paste( dirname(out_file), '/Enrichment/', sep='')
if(!dir.exists(out_dir)){
  dir.create(file.path(out_dir), recursive=T)
}

prey_idx <- grep('^prey$',names(dat.clusters))
bait_idx <- grep('^cluster$',names(dat.clusters))

#res <- Enrichment.main(data_file=dat.clusters, output_dir=out_dir, bait_colname='cluster', prey_colname='prey', grouped=T, enrichment_p_cutoff=1, id_type="UNIPROT")
res <- Enrichment.main(data_file=dat.clusters, output_dir=output_dir, prey_idx=prey_idx, set_idx=set_idx, grouped=T, enrichment_p_cutoff=enrichment_p_cutoff, id_type="UNIPROT", species='mouse')













































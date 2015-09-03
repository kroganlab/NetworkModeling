# Read in a file with a column of metabolite peaks and m/z, search HMDB for similar metabolites based on their weight, adduct, and +/- some threshold.
#   Then compare the associated entrez id's to the metabolites with the significant genes found in the other experiments (RNA seq, PH, UB, etc)
#   Return any significant genes and experiments found to match the metabolites.


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  LOAD FUNCTIONS INTO MEMORY FIRST ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

searchHMDB <- function(x, hmdb, THRESH){
  cat('>> Using fuzzy search to find matching masses in HMDB...\n')
  H= 1.007276
  NH4=18.033823
  Na=22.989218
  K=38.963158
  Cl=34.969402
  
  metabolites <- unique(x$m.z)
  metabolite_matches = list()
  pb <- txtProgressBar(min=0, max=length(metabolites), style=3)
  for(i in 1:length(metabolites)){
    # NOTE: the reported m.z value already includes the adduct, so in order to search the HMDB db, we need to use the straight mass value
    #       without the adduct. EG: mHp is m.z = M+H, so to get M, we need to subtract H from m.z . See below formulas.
    mHp = metabolites[i] - H
    mNH4p = metabolites[i] - NH4
    mNap = metabolites[i] - Na
    mKp = metabolites[i] - K
    mHn = metabolites[i] + H
    mCln = metabolites[i] + Cl
    
    mHp = which( (hmdb$weight >= (mHp-THRESH)) & (hmdb$weight <= (mHp+THRESH)) )
    mNH4p = which( (hmdb$weight >= (mNH4p-THRESH)) & (hmdb$weight <= (mNH4p+THRESH)) )
    mNap = which( (hmdb$weight >= (mNap-THRESH)) & (hmdb$weight <= (mNap+THRESH)) )
    mKp = which( (hmdb$weight >= (mKp-THRESH)) & (hmdb$weight <= (mKp+THRESH)) )
    mHn = which( (hmdb$weight >= (mHn-THRESH)) & (hmdb$weight <= (mHn+THRESH)) )
    mCln = which( (hmdb$weight >= (mCln-THRESH)) & (hmdb$weight <= (mCln+THRESH)) )
    
    idx = c(mNH4p, mNap, mKp, mHn, mCln)
    if( sum(idx)>0 ){
      idx = unique(idx)
      metabolite_matches[[i]] = cbind(hmdb[idx,], "m.z"=metabolites[i])
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  tmp = do.call(rbind,metabolite_matches)
  return(tmp)
}


#  Main function to search handle search
searchDB4HMDB <- function(input_file, hmdb_file, hmdb.entrez_file, flu_file, out_file, THRESH){
  cat(">> Loading files...\n")
  # Load metabolite list
  dat <- input_file
  # Load hmdb
  hmdb <- read.delim(hmdb_file, sep='\t', stringsAsFactors=F, header=T)
  # load hmdb file with corresponding entrez id's
  hmdb.entrez <- read.delim(hmdb.entrez_file, sep='\t', stringsAsFactors=F, header=T)
  names(hmdb.entrez)[grep('entrez_gene_id_mouse',names(hmdb.entrez))] = 'entrez_id'
  
  # Load significant hits from other datasets
  flu <- read.delim(flu_file, sep='\t', stringsAsFactors=F, header=T)
  flu <- flu[,c('experiment_id','omics_type','condition_2','cell_line','strain','entrez_id','symbol')]   # remove unnecessary variables
  
  # search for masses in HMDB 
  hmdb_hits = searchHMDB(dat, hmdb, THRESH)
  cat(">> Matching Peaks to significant hits in other data sets...\n")
  # match up all Peak names with the hits
  hmdb_hits = merge(hmdb_hits, dat[,c('Protein','m.z')], by='m.z')
  # get entrez id's corresponding to HMDB id's
  hmdb_hits = merge(hmdb_hits, hmdb.entrez[,c('hmdb', 'entrez_id')], by.x='id', by.y='hmdb')    ##### for MOUSE
  # check the rest of the flu datasets for significant matches
  hmdb_hits.long = merge(hmdb_hits, flu, by.x='entrez_id', by.y='entrez_id')
  
  # aggregate all the gene names and experiment_id's
  cat(">> Adding gene names and experiment ID's...\n")
  tmp = aggregate(data=hmdb_hits.long[c('id','m.z','Protein','weight','symbol')], symbol~., function(y){ paste(unique(y) ,collapse=",") } )
  tmp2 = aggregate(data=hmdb_hits.long[c('id','m.z','Protein','weight','omics_type')], omics_type~., function(y){ paste(unique(y) ,collapse=",") } )
  # combine the data with the aggregate omics type and symbol
  tmp = merge(tmp, tmp2, by=c('id','m.z','Protein','weight'))
  names(tmp) = c('id','m.z','Protein','matched_weight')
  
  cat(">> Writing results to file...\n")
  write.table(tmp ,gsub(".txt", "_MATCHED.txt", out_file), quote=F, row.names=F, col.names=T, sep='\t')
  write.table(hmdb_hits.long ,gsub(".txt", "_MATCHED_long.txt", out_file), quote=F, row.names=F, col.names=T, sep='\t')
  
  cat(">> Search Complete!!\n")
  return(tmp)
}

#~~~~~~~~~~~~~~~~~~~~~~~~  END FUNCTIONS SECTION  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~










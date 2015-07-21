#! /usr/bin/Rscript
suppressMessages(library(optparse))

mapMasses <- function(results_file, evidence_file, out_file){
  cat('>> MAPPING MASSES\n')
  cat('>>   Reading files...\n')
  dat = read.delim(results_file, stringsAsFactors=F)
  masses = read.delim(evidence_file, stringsAsFactors=F) 
  masses = unique(masses[,c('Proteins','m.z')])
  
  cat('>>   Matching masses...\n')
  ## !!!!!! Keep in mind that there are potetially many masses to 1 sample name due to the rounding used in naming the samples
  mapped = merge(dat, masses, by.x='Protein', by.y='Proteins')
  mapped = mapped[,c(1,length(mapped),3:length(mapped)-1)]
  
  cat('>>   Writing files...\n')
  write.table(mapped, out_file, quote=F, row.names=F, col.names=T, sep="\t")
  return(mapped)
}


annotateMasses <- function(dat,id_file){
  ids = read.delim(id_file, sep='\t', stringsAsFactors=F) 
  tmp = merge ( dat, ids[,-grep('term',names(ids))], by='m.z', all.x=T)
  idx = which(!is.na(tmp$KEGG))
  tmp$Protein[idx] = tmp$KEGG[idx]
  return(tmp[,1:2])
}


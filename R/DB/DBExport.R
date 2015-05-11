source('scripts/DBFunctions.R')

## FIRST OPEN TERMINAL AND TYPE:
## ssh -L 3307:localhost:3306 everschueren@higgs.ucsf.edu

genomics = getGenomics('.*','.*',2,.001, 'genomics_rnaseq')
## doesn't work yet
#genomics_p3 = getGenomics('HTBE','H5N1',2,.001, 'genomics_rnaseq_p3')
proteomics_ph = getModProteomics('.*','.*',1,.05, 'proteomics_ph_sites')
proteomics_ub = getModProteomics('.*','.*',1,.05, 'proteomics_ub_sites')
metabolomics = getMetabolomics('.*','.*',1,.05, 'metabolomics')
basemap = getBaseMap('basemap_menche')

if(!is.null(genomics)){
  write.table(genomics, file='database/tempoffline/genomics.txt', eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}
if(!is.null(proteomics_ph)){
  write.table(proteomics_ph, file='database/tempoffline/proteomics_ph.txt', eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}
if(!is.null(proteomics_ub)){
  write.table(proteomics_ub, file='database/tempoffline/proteomics_ub.txt', eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}
if(!is.null(metabolomics)){
  write.table(metabolomics, file='database/tempoffline/metabolomics.txt', eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}
if(!is.null(basemap)){
  write.table(basemap, file='database/tempoffline/basemap.txt', eol='\n', sep='\t', quote=F, row.names=F, col.names=T)
}
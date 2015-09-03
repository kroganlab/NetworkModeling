
# Source files
source('~/projects/NetworkModeling/R/metabolomics/mapMasses.R')
source('~/projects/NetworkModeling/R/metabolomics/searchDB4HMDB_flat.R')


# ~~~~~~~~~~~~~~~
# MAPPING MASSES
# ~~~~~~~~~~~~~~~
# This section will take an MSStats wide formatted file and map the non-rounded masses to the results for identification purposes.
results_file = '~/Box Sync/Projects/FluomicsModeling//data/datasets-fluomics/metabolomics_colaboration/results/H1N1_lung/H1N1_lung-results-wide.txt'
evidence_file = '~/Box Sync/Projects/FluomicsModeling//data/datasets-fluomics/metabolomics_colaboration/data/H1N1_lung_evidence.txt'
out_file = '~/Box Sync/Projects/FluomicsModeling//data/datasets-fluomics/metabolomics_colaboration/results/H1N1_lung/H1N1_lung-results-wide-masses.txt'
log2FC = 2
pvalue = 0.05

# Mapp the non-rounded masses to the peak names ( rounded m.z./time )
to_search = mapMasses(results_file, evidence_file, out_file)
# Keep only significant fold changes in the results
to_search = to_search[(abs(to_search$LungH1.LungPBS_log2FC)>log2FC) & (to_search$LungH1.LungPBS_adj.pvalue<pvalue),]
# prep for the DB search
to_search = to_search[,1:2]

# ~~~~~~~~~~~~~~~~~~~~~
# ANNOTATE KNOWN PEAKS
# ~~~~~~~~~~~~~~~~~~~~~
{# This section incorporates the id's of any identified metabolites in terms of the KEGG id's so we can use those and not the mz/time
# metabolite_ids = '~/Box Sync/Projects/FluomicsModeling/data/datasets-fluomics/metabolomics_colaboration/Metabolite_IDs/identified_metabolites.txt'
# to_search = annotateMasses(to_search, metabolite_ids)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SEARCH HMDB AGAINST OTHER DATA SETS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This version is designed to work with the flat files offsite (not live mysql db)
# This section will take the non-rounded masses and search for matches in HMDB. The search process includes including adducts and
#   searches for any matches within a window +/- a THRESHold

# database files necessary
hmdb_file = '~/Box Sync/Projects/FluomicsModeling/data/datasets-fluomics/metabolomics_colaboration/Metabolite_IDs/peak_matching/all_hmdb.txt'
flu_file = '~/Box Sync/Projects/FluomicsModeling/data/datasets-fluomics/metabolomics_colaboration/Metabolite_IDs/peak_matching/AllSignificantData.txt'  #includes both human and mouse data, but the ertrez id's won't map between the species so ok to use with both

#!!!!!!!!!! MOSUE or HUMAN??
# select MOUSE or HUMAN db to search
hmdb.entrez_file = '~/Box Sync/Projects/FluomicsModeling/data/datasets-fluomics/metabolomics_colaboration/Metabolite_IDs/peak_matching/HMDB2entrez_mouse.tsv'
#hmdb.entrez_file = '~/Box Sync/Projects/FluomicsModeling/data/datasets-fluomics/metabolomics_colaboration/Metabolite_IDs/peak_matching/HMDB2entrez_human.tsv'

THRESH = 0.05  # this is the amount we are willing to let the masses be off for identification +/-

# Search the Fluomics database for genes associated with metabolites
#~~~~~~~~~~~~~~~~~~~~
results = searchDB4HMDB(to_search, hmdb_file, hmdb.entrez_file, flu_file, out_file, THRESH)




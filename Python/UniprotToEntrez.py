#!/usr/bin/env python
 
import pandas
import mygene
 
mg = mygene.MyGeneInfo()
 
# protein_table = pandas.read_table('/Users/mchang/Downloads/human_uniprot_20150115.txt', index_col=0)
protein_table = pandas.read_table('/Users/everschueren/Projects/HPCKrogan/Data/FluOMICS/database/resources/human_uniprot_20150115.txt', index_col=0)
 
q = protein_table.index.tolist()
# q.append('Q5JQC4')
 
# uniprot_entrez = mg.querymany(q, scopes="uniprot", species="human", fields="entrezgene", as_dataframe=True, df_index=False, verbose=False)
uniprot_entrez = mg.querymany(q, scopes="uniprot", species="human", fields="entrezgene", as_dataframe=True, df_index=False, verbose=False)
 
# get rid of rows without Entrez Gene ID and convert to string
uniprot_entrez.dropna(inplace=True, subset=["entrezgene"])
uniprot_entrez['entrezgene'] = uniprot_entrez['entrezgene'].astype(int)
uniprot_entrez['entrezgene'] = uniprot_entrez['entrezgene'].apply(str)
 
# build semicolon-delimited Gene IDs by uniprot ID (query)
uniprot_group_entrez = uniprot_entrez.groupby('query')
uniprot_entrez_series = uniprot_group_entrez.entrezgene.apply(lambda x: ';'.join(sorted(x)))
 
# merge new Gene IDs into table and write out
add = protein_table.join(uniprot_entrez_series)
# add.to_csv("human_uniprot_test.csv")
add.to_csv("human_uniprot_.csv")
 
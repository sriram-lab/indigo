"""
Pseudocode

1. Remove all special characters from genes of new list
2. Remove beginning portion of gene names from INDIGO chemogenomics
2. Figure out which genenames match, and rearrange new gene list with all data in order of ecoli phenotype file
3. Put all remaining generows at the end. Might want it to have same structure as ecoli, maybe need to use David to look up genes
and get that first identiifer part.
"""

#could always use ismember
# want to get genes in new data that are in old data and where they are located in the old data so that you can rearrange them in that way
# and then whatever elements are not there, just append them to the end.
# prob much easier in matlab

import pandas as pd

path = '/Users/david/Documents/systems-bio-research/indigo/data/maeda_files_ecoli_transcriptomics/'

# New gene list
maeda_data = pd.read_excel(path + 'maeda_processed.xlsx', sheet_name='maeda_log2FC')

maeda_genes = maeda_data['gene']

# process gene list
# remove all non alpha numeric values

# Original gene list

orig_chemgen = pd.read_excel(path+'orig_chemgen.xlsx')
orig_genes = orig_chemgen['Gene']
print(orig_genes)
# print(maeda_genes.intersection(orig_genes))
import pandas as pd

# only include drugs that are in good drug list in final output
data = pd.read_excel('pairwise_interactions_processed.xlsx',header=None)
good_drugs = ['CHL','CYS','FUS','INH','NIT','RIF','TUN','VAN']

data = data[data.iloc[:,0].isin(good_drugs)]
data = data[data.iloc[:,1].isin(good_drugs)]
data.to_excel('pairwise_interactions_final.xlsx')
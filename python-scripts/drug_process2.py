import pandas as pd

# only include drugs that are in good drug list in final output
data = pd.read_excel('mtb_yilancioglu_2019_complete.xlsx',header=None)

#this can be modified
good_drugs = ['CHL','CYS','FUS','INH','NIT','RIF','TUN','VAN']
# good_drugs = ['BDQ','CHL','CLZ','CYS','ETH','ETA','FUS','INH','LIN','MOX','NIT','RIF','SQ1','TUN','VAN']

data = data[data.iloc[:,0].isin(good_drugs)]
data = data[data.iloc[:,1].isin(good_drugs)]
data.to_excel('pairwise_interactions_final.xlsx',header=False, index=False)
import pandas as pd

data = pd.read_excel('pairwise_interactions.xlsx')

# iterate through each column
# iterate through the rows for each column
# if there is a 1, assign it the value of the column name
for i in data.columns:
    data[i].replace(1,i.upper(),inplace=True)

data.to_excel('pairwise_interactions_processed.xlsx')
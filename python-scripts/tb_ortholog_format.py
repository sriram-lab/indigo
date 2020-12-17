import pandas as pd
import os

dir = '/Users/david/Downloads/a'

for file in os.listdir(dir):
    if file.endswith('.xlsx'):
        file = dir + '/' + file
        df = pd.read_excel(file,header=None)
        df.iloc[:,0] = df.iloc[:,0].str.replace('BD_','')
        df.iloc[:,0] = df.iloc[:,0].str.slice_replace(start=0,stop=2,repl='Rv')
        with pd.ExcelWriter(file) as writer:  
            df.to_excel(writer,index=False,header=False)
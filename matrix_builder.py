"""The script builds COG matrix from database"""
import pandas as pd

cogdb_path = 'data/cog2003-2014.csv'
cog_matrix_output_path = 'data/cog_matrix.csv'

db = pd.read_csv(cogdb_path, header=None, skipinitialspace=True)
data = {}
for genome in db[1].unique():
    data[genome] = {}
    for cog in db[db[1] == genome][6].unique():
        data[genome][cog] = 1
print('Done!')
pd.DataFrame.from_dict(data, dtype=int).fillna(0).to_csv(cog_matrix_output_path)

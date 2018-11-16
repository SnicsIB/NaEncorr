import pandas as pd
from scipy.spatial.distance import pdist, squareform

be_path = 'data/bioenerg.csv'
matrix_path = 'data/cog_matrix.csv'
scope_path = 'data/membrane.txt'  # list of cogs for investigation
names_path = 'data/cognames2003-2014.tab'
corr_matrix_output = 'data/corrs.csv'

# bioenergetics loading
bioenergetics = pd.read_csv(be_path, index_col=6)
bioenergetics['Coupling ion'].replace(
    to_replace=['Na+', 'H+', 'Na+ (H+, N)', 'Na+ (H+)', 'H+ (N)'],
    value=[1, 0, 0, 0, 0],
    inplace=True
)

# input scope
scope = list(map(lambda x: x.strip(), open(scope_path)))

# input names
names = pd.read_table(names_path, index_col=0)

# input of matrix and filtering
df = pd.read_csv(matrix_path, index_col=0)
df = df.loc[scope]

# adding description to COGs labels
labels = ['{}; {}'.format(x, names.loc[x]['name']) for x in scope]
df.index = labels

df.loc['Sodium bioenergetics'] = bioenergetics['Coupling ion']
df.dropna(inplace=True, axis=1)
# df = df[[x for x in df.columns if bioenergetics.loc[x]['Kingdom'] == 'Archaea']] # you can change scope HERE!
df = df[(df.T != 0).any()]  # delete empty rows
df = df[df.T.mean() < 0.6]  # make it sparse (too frequent COGs' removal)

# computing correlations
corrs = pd.DataFrame(1 - squareform(pdist(df.values, metric='dice')))
corrs.columns = list(df.index)
corrs.index = list(df.index)
corrs.sort_values('Sodium bioenergetics').to_csv(corr_matrix_output)

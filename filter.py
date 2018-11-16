"""Script for choosing only membrane cogs"""
import seaborn as sns
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt

tmhmm_path = 'data/711Genomes.tmhmm.txt'
cogdb_path = 'data/cog2003-2014.csv'
threshold = 0.99
profile_img_output = 'profile.svg'
membr_cogs_list_output = 'data/membrane.txt'


# proteins with two or more predicted helices
membrane_prots = [
    t.split('|')[0] for t in list(map(lambda x: x.strip(), open(tmhmm_path)))
    if int(t.split('PredHel=')[1].split('\t')[0]) > 1
]

db = pd.read_csv(cogdb_path, header=None, skipinitialspace=True)
membr_prots_in_cog = Counter(db[db[0].isin(membrane_prots)][6])
prots_in_cog = Counter(db[6])

membrane_part = {}
for key in membr_prots_in_cog:
    membrane_part[key] = membr_prots_in_cog[key]/prots_in_cog[key]

# drawing distribution of membranous part in cogs
plt.figure(figsize=(12, 9))
plt.title('Membrane proteins parts in COGs. Distribution')
plt.xlabel('membrane proteins part')
sns.distplot(list(membrane_part.values()), kde=None, bins=100)
plt.tight_layout()
plt.savefig(profile_img_output)

# writing COGs with MP >= threshold
with open(membr_cogs_list_output, 'w') as fout:
    for key in membrane_part:
        if membrane_part[key] >= threshold:
            fout.write(key+'\n')


# coding: utf-8

# In[47]:


import os
from iteration_utilities import deepflatten, Iterable
import pandas as pd
from Bio import SeqIO

wdir = "711_genomes"
interesting_cogs = set([
    "COG1822",
    "COG1906",
    "COG2456",
    "COG2426",
    "COG4035",
    "COG4720",
    "COG4769",
    "COG2245",
    "COG1967"
])


def get_location(description):
    """be careful: cannot resolve around 1st nucleotide"""
    res = dict()
    res["strand"], res["loc"] = description.split("|")[4:6]
    res["strand"] = int(res["strand"])
    loc = res["loc"].split("..")
    try:
        l1 = int(loc[0])
    except ValueError:
        l1 = int(loc[0][1:])
    try:
        l2 = int(loc[-1])
    except ValueError:
        l2 = int(loc[-1][1:])
    res["loc"] = (l1, l2)
    return res
    
        


def get_operons(seqio_handler, max_dist=80):
    """accepts sorted by coordinates genes in fasta format. 5 - strand, 6 - coordinates"""
    operons_list = [[]]
    current_operon = -1
    old_location = location = {"strand": 0, "loc":(0, 0)}
    for seq in seqio_handler:
        location = get_location(seq.description)
        if location["strand"] == old_location["strand"] and location["loc"][0] - old_location["loc"][1] < max_dist:
            operons_list[current_operon].append(seq)
            # if location["loc"][0] - old_location["loc"][1] < 0:  # check bypassing genome's first nucleotide
                # print("Warning! 1 nucleotide passed.")
        else:
            current_operon += 1
            operons_list.append([])
            operons_list[current_operon].append(seq)
        old_location = location
    return operons_list


def gi_to_cog(df, gi):
    try:
        series = df.loc[gi][1]
        if isinstance(series, str):
            return [series]
        else:
            return list(series)
    except KeyError:
        return 0


# In[48]:


res_dict = {}
total_operons = {}
all_closest_neighbours = {} # only for insteresting cog
for folder in filter(lambda x: os.path.isdir(wdir+"/"+x), os.listdir(wdir)):
    print("Processing {}...".format(folder))
    operons = get_operons(SeqIO.parse("{}/{}/p_{}.fasta".format(wdir, folder, folder), "fasta"))  # get grouped sequences
    try:
        gi_to_cogs = pd.read_table("{}/{}/c_{}.txt".format(wdir, folder, folder), index_col=0, header=None) # table of gis
    except:
        gi_to_cogs = pd.DataFrame({0:[1]})  # if table is empty
    # Seq to COGs
    cog_operons = []
    for t in operons:
        cog_operons.append([])
        for x in t:
            xcogs = gi_to_cog(gi_to_cogs, int(x.id.split("|")[1]))
            if xcogs != 0:
                for xcog in xcogs:
                    if xcog in interesting_cogs:
                        if xcog not in all_closest_neighbours.keys():
                            all_closest_neighbours[xcog] = t[:]
                        else:
                            all_closest_neighbours[xcog] += t
            cog_operons[-1].append(xcog)
    # constructing table as dict (count collocations of every pair of COGs)
    for operon in cog_operons: 
        for pair in Iterable(deepflatten(operon, ignore=str)).combinations(2):
            if pair[0] not in res_dict.keys():
                res_dict[pair[0]] = {}
            if pair[1] not in res_dict[pair[0]].keys():
                res_dict[pair[0]][pair[1]] = 0
            res_dict[pair[0]][pair[1]] += 1
            # symmetric pair
            if pair[1] not in res_dict.keys():
                res_dict[pair[1]] = {}
            if pair[0] not in res_dict[pair[1]].keys():
                res_dict[pair[1]][pair[0]] = 0
            res_dict[pair[1]][pair[0]] += 1
            # addd one operon to counter
            if pair[1] not in total_operons.keys():
                total_operons[pair[1]] = 0
            if pair[0] not in total_operons.keys():
                total_operons[pair[0]] = 0
            total_operons[pair[1]] += 1
            total_operons[pair[0]] += 1
# [index cog] needs [column cog] coefficient -> 1
for icog in res_dict:
    for ccog in res_dict[icog]:
        res_dict[icog][ccog] /= total_operons[icog]


# In[49]:


for key in all_closest_neighbours:
    SeqIO.write(all_closest_neighbours[key], "{}.fasta".format(key), "fasta")
res_df = pd.DataFrame(res_dict)
res_df.fillna(0, inplace=True)
res_df.to_csv("conditional.csv") # output raw data


# In[50]:


# name cogs
cog_names = pd.read_table("cognames2003-2014.tab", index_col=0)
res_dft = res_df.join(cog_names["name"])
# print best results for chosen COGs
for cog in interesting_cogs:
    try:
        res_dft.sort_values(cog, ascending=False).head(20)[["name", cog]].to_csv("besthits_{}.csv".format(cog))
    except KeyError:
        print("No hits for {}".format(cog))
# print best results for chosen COGs (reversed hits)
res_dft = res_df.T
res_dft = res_dft.join(cog_names["name"])
for cog in interesting_cogs:
    try:
        res_dft.sort_values(cog, ascending=False).head(20)[["name", cog]].to_csv("rbesthits_{}.csv".format(cog))
    except KeyError:
        print("No hits for {}".format(cog))

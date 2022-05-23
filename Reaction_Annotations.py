
import re
import requests
import pandas as pd
from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser


all_reactions = pd.read_csv("reactionpresence.csv", nrows=1).columns.tolist()
all_reactions_strip = [x[:-2] for x in all_reactions]
all_reaction_dict = dict.fromkeys(all_reactions_strip)



r = requests.get("https://raw.githubusercontent.com/ModelSEED/ModelSEEDDatabase/master/Biochemistry/reactions.tsv")
r = r.content.decode("utf-8")

count = 0
for key in all_reaction_dict:
    annotations = r
    present = annotations.find(key)
    if present != -1:
        right = annotations[present:]
        position = right.find("KEGG: ")
        if position != -1:
            all_reaction_dict[key] = right[position+6:position+13]
        elif position == -1:
            all_reaction_dict[key] = "NA"
    elif present == -1:
        all_reaction_dict[key] = "NA"
    print(count)
    count+=1
clean_dict = {key.strip("|"): item.strip("|") for key, item in all_reaction_dict.items()}


count = 0
for key, value in clean_dict.items():
    try:
        kegg = REST.kegg_get(value).read()
        position = kegg.find("Metabolism; ")
        if position != -1:
            clean_dict[key] = re.search('Metabolism;(.*)\n', kegg)
        elif position == -1:
            clean_dict[key] = "NA"
    except:
        pass
    print(count)
    count +=1


with open('reaction_annotations.csv', 'w') as f:
    count = 0
    for key in clean_dict.keys():
        f.write("%s,%s\n"%(key,clean_dict[key]))
        print(count)
        count +=1







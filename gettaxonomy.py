import csv
from ete3 import NCBITaxa
import pandas as pd
import cobra
import glob
import natsort

my_file = open("taxids.txt", "r")
data = my_file.read()
taxids = data.split("\n")

my_file = open('genomeids.txt', "r")
data = my_file.read()
patricids = data.split('\n')

ncbi = NCBITaxa()

def get_desired_ranks(taxid, desired_ranks):
    lineage = ncbi.get_lineage(taxid)   
    names = ncbi.get_taxid_translator(lineage)
    lineage2ranks = ncbi.get_rank(names)
    ranks2lineage = dict((rank,taxid) for (taxid, rank) in lineage2ranks.items())
    return{'{}_id'.format(rank): ranks2lineage.get(rank, '<not present>') for rank in desired_ranks}

def get_model_info():
    filelist = glob.glob("*.sbml")
    filelist = natsort.natsorted(filelist)
    numrxns = []
    numgenes =[]
    nummetabol = []
    strain = 0
    for file in filelist:
        model = cobra.io.read_sbml_model(file, low_memory=False)
        numrxns.append(len(model.reactions))
        numgenes.append(len(model.genes))
        nummetabol.append(len(model.metabolites))
        strain += 1
        print(strain)
    return numrxns, numgenes, nummetabol



if __name__ == '__main__':
    desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    results = list()
    for taxid in taxids:
        results.append(list())
        results[-1].append(str(taxid))
        ranks = get_desired_ranks(taxid, desired_ranks)
        for key, rank in ranks.items():
            if rank != '<not present>':
                results[-1].append(list(ncbi.get_taxid_translator([rank]).values())[0])
            else:
                results[-1].append(rank)
    #get model info
    numrxns, numgenes, nummetabol = get_model_info()
    
    #generate the header
    header = ['Original_query_taxid']
    header.extend(desired_ranks)
    print('\t'.join(header))

    #print the results
    resultarr= []
    for result in results:
        print('\t'.join(result))
        resultarr.append(result[1])

    taxonomyinfo = pd.DataFrame({'PATRIC ID': patricids, 'Taxid': taxids, 'Phylum': resultarr, 'Reactions': numrxns, 'Genes': numgenes, 'Metabolites': nummetabol})
    taxonomyinfo.to_csv('taxonomyinfo.csv')




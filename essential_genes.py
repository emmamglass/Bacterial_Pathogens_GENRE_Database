import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cobra
from cobra.io import read_sbml_model
import glob, os
import natsort
import ast
import json

#Determine universal essential genes across models
#create empty dictionary, keys are PATRIC IDs

def find_all_essential_genes():
	keys = glob.glob('*.sbml')
	keys = natsort.natsorted(keys)

	mother_dict = dict.fromkeys(keys)
	print(mother_dict)

	#determine core essential genes per strain, save set of essential genes
	#to the value associated with PATRIC ID key
	i = 0
	for key in mother_dict:
		i +=1
		print(i)
		file = str(key)
		model = read_sbml_model(file)
		essential_genes = cobra.flux_analysis.variability.find_essential_genes(model)
		mother_dict[key] = essential_genes
		print(mother_dict)

	return mother_dict

def find_core_essential_genes():



	core = reduce(set.intersection, (set(val) for val in mother_dict.values()))

	print(core)

	return core 

def read_essential_genes():

	with open("allessentialgenes.txt","r") as data:
		data = str(data.read())
		data = data.replace("'", '"')
		print(data)

		dictionary = json.loads(data)
	return 


'''def find_core_essential_genes_by_class(mother_dict):
	patric_ids = pd.read_excel('Taxids.xlsx')
	classes = pd.read_excel('classchron.xlsx')

	id_class = pd.concat(patric_ids, classes)'''






 

if __name__ == '__main__':
	read_essential_genes()





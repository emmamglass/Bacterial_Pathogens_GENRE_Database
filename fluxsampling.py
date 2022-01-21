import cobra.test
from cobra.io import read_sbml_model
import os
import glob
from os.path import join
from cobra.sampling import sample
import pandas as pd
import numpy as np
from numpy.linalg import svd
from sklearn.decomposition import TruncatedSVD
from matplotlib import pyplot as plt
from sklearn import manifold
from sklearn.metrics import euclidean_distances
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist
from sklearn.metrics.pairwise import pairwise_distances
import argparse


parser = argparse.ArgumentParser(description = 'Reads .sbml models, performs flux sampling, and plots flux data, clustering on a certian input (pathogen name, oxygen usage status, gram')
parser.add_argument('--cluster_type', default = 'none')
parser.add_argument('--sample_size', default = 'none')
args = parser.parse_args()

#Flux sampling functions
def get_bacteria_names():
	sbml_files = glob.glob(('*.sbml'))
	bacteria = []
	for file in sbml_files:
		str(file)
		bac_name = file.split('.')
		bac_name.pop()
		bac_name.pop()
		bacteria += bac_name
	return bacteria

def read_recons():
	#gets the names of all bacteria that there are sbml reconstructions avaliable for
	#creates a list of all bacteria names
	#get_bacteria_names()
	sbml_files = glob.glob(('*.sbml'))
	bacteria = []
	for file in sbml_files:
		str(file)
		bac_name = file.split('.')
		bac_name.pop()
		bac_name.pop()
		bacteria += bac_name
	
	#opens all sbml models into a list 
	models = []
	for file in glob.glob("*.sbml"):
		models.append(read_sbml_model(file, low_memory=False))

	#create dictionary of reconstructions where the keys are the bacteria names and the 
	#values are the sbml reconstructions
	recons = {k:v for k, v in zip(bacteria, models)}
	return recons

def sample_flux(recons, sample_size):
	#for every sbml model, flux is sampled for the specified sample size
	for bacteria, model in recons.items():
		flux = sample(model, sample_size)
		flux.to_csv(str(bacteria) + ".csv")
	return 




#NMDS functions
def read_flux_data():
	#get all bacteria names from csv files
	flux_csv = glob.glob(('*.csv'))
	bacteria = []
	for file in flux_csv:
		str(file)
		bac_name = file.split('.')
		bac_name.pop()
		bacteria += bac_name

	#read in all csv files into a list
	data = []
	for file in glob.glob('*.csv'):
		data.append(pd.read_csv(file))

	#create a dictionary with bacteria names and flux data
	flux = {k:v for k, v in zip(bacteria, data)}
	return flux

def add_pathogen_data(flux):
	pathogen_data = []
	for bacteria, flux_data in flux.items():
		flux_data['pathogen'] = str(bacteria).replace("_",' ')
		pathogen_data.append(flux_data)
	return pathogen_data

def add_gram_data(flux):
	gram_data = []
	for bacteria, flux_data in flux.items():
		if bacteria not in ('Abiotrophia_defectiva', 'Acidaminococcus_fermentans', 'Aerococcus_viridians', 'Amycolatopsis_orientalis', 'Eubacterium_brachy'):
			flux_data['gram'] = 'gram negative'
		else:
			flux_data['gram'] = 'gram positive'
		gram_data.append(flux_data)
	return gram_data

def add_aerobic_data(flux):
	aerobic_data = []
	for bacteria, flux_data in flux.items():
		if bacteria not in ('Abiotrophia_defectiva', 'Actinomyces_gerencseriae', 'Aerococcus_viridians','Acidaminococcus_fermentans', 'Bacteroides_uniformis', 'Eubacterium_brachy'):
			flux_data['aerobic'] = 'aerobe'
		elif bacteria not in  ('Abiotrophia_defectiva', 'Actinomyces_gerencseriae', 'Aerococcus_viridians', 'Amycolatopsis_orientalis', 'Acinetobacter_baumannii', 'Achromobacter_xylosoxidans', 'Achromobacter_piechaudii') :
			flux_data['aerobic'] = 'anaerobe'
		else:
			flux_data['aerobic'] = 'facultative anaerobe'
		aerobic_data.append(flux_data)
	return aerobic_data


def consensus_rxns(data):
	#print(data[0].columns)
	#returns a large data frame of all consensus reactions for all samples of all bacteria
	rxns = []
	for i in range(len(data)-1):
		rxns += [col for col in set(data[i].columns).intersection(data[i+1].columns)]

	rxn_list = []
	count_list = []
	for rxn in rxns:
		rxn_list.append(rxn)
		counts = rxns.count(rxn)
		count_list.append(counts)

	rxn_counts = {k:v for k, v in zip(rxn_list, count_list)}

	consensus_rxns = []
	for rxn, count in rxn_counts.items():
		if count == len(data)-1:
			consensus_rxns.append(rxn)

	processed_data = []
	for bacteria in data:
		new_data = bacteria[consensus_rxns]
		processed_data.append(new_data)

	all_data = pd.concat(processed_data, ignore_index = True)

	if 'pathogen' in all_data:
		cluster_data = all_data[['pathogen']]
		del all_data['pathogen']
	elif 'gram' in all_data:
		cluster_data = all_data[['gram']]
		del all_data['gram']
	elif 'aerobic' in all_data:
		cluster_data =all_data[['aerobic']]
		del all_data['aerobic']

	return all_data, cluster_data


def perform_NMDS(data):
	nmds = MDS(n_components = 2, n_init =1, metric = True, random_state = 0, dissimilarity = 'precomputed')

	similarities = pairwise_distances(data, metric = 'braycurtis')
	data_transformed = nmds.fit_transform(similarities)

	return data_transformed

def plot_pathogen(data, pathogen_data, bacteria):
	df = pd.DataFrame(data)
	df = pd.concat([pathogen_data, df], axis = 1)
	df = df.apply(lambda x: pd.Series(x.dropna().values)).fillna('')
	
	targets = ['Abiotrophia defectiva', 'Achromobacter piechaudii', 'Achromobacter xylosoxidans', 'Acinetobacter baumannii', 'Acidaminococcus fermentans', 'Actinomyces gerencseriae', 'Aerococcus viridians', 'Amycolatopsis orientalis', 'Bacteroides uniformis', 'Eubacterium brachy']
	colors = ['lightcoral', 'sandybrown','gold','yellowgreen', 'mediumseagreen', 'mediumturquoise', 'cornflowerblue', 'darkorchid', 'orchid', 'slategrey']

	for target, color in zip(targets, colors):
		indicesToKeep = df['pathogen'] == target
		plt.scatter(df.loc[indicesToKeep,0],df.loc[indicesToKeep,1], c = color, alpha = 0.7)

	plt.legend(bacteria)

	plt.xlabel('NMDS1')
	plt.ylabel('NMDS2')

	plt.title('NMDS')

	plt.show()

	return

def plot_gram(data, gram_data):
	df = pd.DataFrame(data)
	df = pd.concat([gram_data, df], axis = 1)
	df = df.apply(lambda x: pd.Series(x.dropna().values)).fillna('')

	targets = ['gram negative', 'gram positive']
	colors = ['dodgerblue', 'indigo']

	for target,color in zip(targets,colors):
		indicesToKeep = df['gram'] == target
		plt.scatter(df.loc[indicesToKeep,0],df.loc[indicesToKeep,1],c=color, alpha = 0.7)

	plt.legend(['gram negative', 'gram positive'])

	plt.xlabel('NMDS1')
	plt.ylabel('NMDS2')

	plt.title('NMDS')

	plt.show()

	return

def plot_aerobic(data, aerobic_data):
	df = pd.DataFrame(data)
	df = pd.concat([aerobic_data, df], axis = 1)
	df = df.apply(lambda x: pd.Series(x.dropna().values)).fillna('')

	targets = ['aerobe', 'anaerobe', 'facultative anaerobe']
	colors = ['darkolivegreen', 'goldenrod', 'turquoise']

	for target,color in zip(targets,colors):
		indicesToKeep = df['aerobic'] == target
		plt.scatter(df.loc[indicesToKeep,0],df.loc[indicesToKeep,1],c=color, alpha = 0.7)

	plt.legend(['gram negative', 'gram positive'])

	plt.xlabel('NMDS1')
	plt.ylabel('NMDS2')

	plt.title('NMDS')

	plt.show()

	return


#################

#process input settings
cluster_type = str(args.cluster_type)
sample_size = int(args.sample_size)
#bypass = str(args.bypass)

bacteria = get_bacteria_names()
recons = read_recons()
flux_sample_data = sample_flux(recons, sample_size)
flux_data = read_flux_data()

if cluster_type == 'gram':
	data = add_gram_data(flux_data)
	processed_data, cluster_data = consensus_rxns(data)
	NMDS_data = perform_NMDS(processed_data)
	plot_gram(NMDS_data, cluster_data)

elif cluster_type == 'pathogen':
	data = add_pathogen_data(flux_data)
	processed_data, cluster_data = consensus_rxns(data)
	NMDS_data = perform_NMDS(processed_data)
	plot_pathogen(NMDS_data, cluster_data, bacteria)

elif cluster_type == 'aerobic':
	data = add_aerobic_data(flux_data)
	processed_data, cluster_data = consensus_rxns(data)
	NMDS_data = perform_NMDS(processed_data)
	plot_aerobic(NMDS_data, cluster_data)













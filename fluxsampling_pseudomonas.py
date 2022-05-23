#import cobra.test
from cobra.io import read_sbml_model
from skbio.stats.distance import permanova
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
from skbio.stats.distance import DistanceMatrix
from sklearn.cluster import KMeans
from skbio.stats.distance import anosim
import argparse



#Parse the users input arguments
parser = argparse.ArgumentParser(description = 'Reads .sbml models, performs flux sampling, and plots flux data, clustering on a certian input (pathogen name, oxygen usage status, gram')
parser.add_argument('--cluster_type', default = 'none')
parser.add_argument('--sample_size', default = '500')
parser.add_argument('--bypass', default = 'false')
parser.add_argument('--consensus_rxns', default = 'false')
args = parser.parse_args()


#Flux sampling functions

#get names of bacterial SPECIES from the .sbml file names
def get_bacteria_names():
	sbml_files = glob.glob(('*.sbml'))
	bacteria = []
	for file in sbml_files:
		str(file)
		bac_name = file.split('.')
		bac_name.pop()
		bacteria += bac_name
		print(bacteria)
	return bacteria


#read all .sbml files in a folder
def read_recons():
	#gets the names of all bacteria that there are sbml reconstructions avaliable for
	#creates a list of the bacteria names
	sbml_files = glob.glob(('*.sbml'))
	bacteria = []
	for file in sbml_files:
		str(file)
		bac_name = file.split('.')
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
	#for every sbml model, flux is sampled for the specified sample size inputted by user
	#flux is saved to 
	for bacteria, model in recons.items():
		#print(bacteria)
		print(bacteriaflux)
		flux = sample(model, sample_size)
		flux.to_csv(str(bacteria) + ".csv")
	return 


#NMDS functions

#reads the flux data from csv file
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

#add a column to the flux data output that specifies the pathogen names
def add_pathogen_data(flux):
	pathogen_data = []
	for bacteria, flux_data in flux.items():
		flux_data['pathogen'] = str(bacteria).replace("_",' ')
		pathogen_data.append(flux_data)
	return pathogen_data


#determine what reactions are present in ALL reconstructions
def consensus_rxns(data):
	#determine what rxns are present in all reconstructions

	#determine consensus reactions in adjacent reconstructions
	rxns = []
	for i in range(len(data)-1):
		rxns += [col for col in set(data[i].columns).intersection(data[i+1].columns)]

	#add adjacent reactions to a list, count how many times a given rxn is in the list
	rxn_list = []
	count_list = []
	for rxn in rxns:
		rxn_list.append(rxn)
		counts = rxns.count(rxn)
		count_list.append(counts)

	#add the number of reaction counts and name of the reaction to a dictionary
	rxn_counts = {k:v for k, v in zip(rxn_list, count_list)}

	#if the reaction appears the same number of times that there are reconstructions
	#in the file, add that reaction to a list of consensus reactions
	consensus_rxns = []
	for rxn, count in rxn_counts.items():
		if count == len(data)-1:
			consensus_rxns.append(rxn)

	#for every reconstruction, keep only the columns that correspond to consensus reactions
	processed_data = []
	for bacteria in data:
		new_data = bacteria[consensus_rxns]
		processed_data.append(new_data)

	#concatenate all dataframes to be able to run into the NMDS function
	all_data = pd.concat(processed_data, ignore_index = True)

	#delete the column specifying the cluster method, but save it to a variable
	if 'pathogen' in all_data:
		cluster_data = all_data[['pathogen']]
		del all_data['pathogen']
	elif 'gram' in all_data:
		cluster_data = all_data[['gram']]
		del all_data['gram']
	elif 'aerobic' in all_data:
		cluster_data =all_data[['aerobic']]
		del all_data['aerobic']

	#return all data used for NMDS, and the vector cluster data that will 
	#determine how we cluster our data
	return all_data, cluster_data

def all_rxns(data):
	#print(data)
	#get rxns for all 
	rxns = []
	for i in range(len(data)-1):
		rxns += [col for col in set(data[i].columns)]

	#count how many times a given rxn is in the list
	rxn_list = []
	count_list = []
	for rxn in rxns:
		rxn_list.append(rxn)
		counts = rxns.count(rxn)
		count_list.append(counts)

	#add the number of reaction counts and name of the reaction to a dictionary
	rxn_counts = {k:v for k, v in zip(rxn_list, count_list)}

	#
	#print(data)
	for rxn, counts in rxn_counts.items():
		for i in range(len(data)-1):
			#print(data[i])
			if str(rxn) not in data[i].columns:
				data[i][str(rxn)] = 0
	
	all_data = pd.concat(data, ignore_index = True).fillna(0)

	if 'pathogen' in all_data:
		cluster_data = all_data[['pathogen']]
		del all_data['pathogen']
	elif 'gram' in all_data:
		cluster_data = all_data[['gram']]
		del all_data['gram']
	elif 'aerobic' in all_data:
		cluster_data =all_data[['aerobic']]
		del all_data['aerobic']

	return all_data, cluster_data, 

#this performs non metric multidimensional scaling
def perform_NMDS(data, cluster_data):
	nmds = MDS(n_components = 2, n_init =1, metric = True, random_state = 0, dissimilarity = 'precomputed')

	similarities = pairwise_distances(data, metric = 'braycurtis')

	'''print(np.shape(similarities))
	cluster_data = cluster_data.to_numpy()
	for i in range(len(cluster_data)):
		if cluster_data[i] == '287':
			cluster_data[i] = 'Non-CF'
		if cluster_data[i] == '1407059':
			cluster_data[i] = 'Non-CF'
		if cluster_data[i] == '1089456':
			cluster_data[i] = 'Non-CF'
		if cluster_data[i] == '1427342':
			cluster_data[i] = 'CF'
		if cluster_data[i] == '1437873':
			cluster_data[i] = 'Non-CF'
		if cluster_data[i] == '381754':
			cluster_data[i] = 'Non-CF'
		if cluster_data[i] == '910265':
			cluster_data[i] = 'CF'
		if cluster_data[i] == '1340851':
			cluster_data[i] = 'CF'
		if cluster_data[i] == '208963':
			cluster_data[i] = 'Non-CF'
		if cluster_data[i] == '1280938':
			cluster_data[i] = 'Non-CF' 

	print(cluster_data)

	permanova_result = permanova(distance_matrix = DistanceMatrix(similarities), grouping = cluster_data) #['910265', '381754', '287', '208963', '1437873', '1427342', '1407059', '1352355', '1340851', '1280938', '1089456'],)
	print(permanova_result)

	anosim_result = anosim(distance_matrix = DistanceMatrix(similarities), grouping = cluster_data)
	print(anosim_result)'''

	data_transformed = nmds.fit_transform(similarities)

	return data_transformed

#This plots the results of the NMDS and clusters based upon pathogen name
def plot_pathogen(data, pathogen_data, bacteria):
	df = pd.DataFrame(data)
	df = pd.concat([pathogen_data, df], axis = 1)
	df = df.apply(lambda x: pd.Series(x.dropna().values)).fillna('')

	print(df['pathogen'].to_string())
	targets = ['1440052', '216592', '331112', '168807', '941323', '419947', '1346787', '1773', '1806', '1352598', '272944', '35790', '1290428', '392021', '1003201', '287', '1340851', '381754', '1407059', '910265']
	#targets = ['Non-CF', 'CF']
	#legend_names = ['Non-CF', 'CF']
	#colors = ['#18851a', '#14a395']
	legend_names = ['P. a. 1', 'P. a. 2', 'P. a. 3', 'P. a. 4', 'P. a. 5', 'M. t. 1', 'M. t. 2', 'M. t. 3', 'M. t. 3', 'M. t. 4', 'M. t. 5', 'R. r. 1', 'R. r. 2', 'R. r. 3', 'R. r. 4', 'R. r. 5', 'P. a. 1', 'P. a. 2', 'P. a. 3', 'P. a. 4', 'P. a. 5']
	colors = ['#75670d', '#a18f1d', '#bdab3e', '#d4c35b','#eddd7e', '#47b562','#32a852', '#248a3f', '#196b2e', '#0f4f1f', '#0f4f4c', '#1b6965', '#26827d', '#3ca39e', '#4cc2bc', '#102b4f', '#2d5891', '#4076bd', '#5089d4', '#78aef5']

	for target, color in zip(targets, colors):
		indicesToKeep = df['pathogen'] == target
		plt.scatter(df.loc[indicesToKeep,0],df.loc[indicesToKeep,1], c = color, alpha = 0.7)

	NMDS1 = df.loc[indicesToKeep,0]
	NMDS2 = df.loc[indicesToKeep,1]

	plt.legend(legend_names)

	plt.xlabel('NMDS1')
	plt.ylabel('NMDS2')

	#plt.title('NMDS')

	#plt.savefig('pathogen.png')
	plt.show()
	return



#
#####################################################################################


#process input settings
cluster_type = str(args.cluster_type)
sample_size = int(args.sample_size)
bypass = str(args.bypass)
consensus_reaxns = str(args.consensus_rxns)

#Run main

if bypass == 'true':
	bacteria = get_bacteria_names()
	flux_data = read_flux_data()
else:
	bacteria = get_bacteria_names()
	recons = read_recons()
	flux_sample_data = sample_flux(recons, sample_size)
	flux_data = read_flux_data()

if cluster_type == 'pathogen':
	data = add_pathogen_data(flux_data)
	if consensus_reaxns == 'true':
		processed_data, cluster_data = consensus_rxns(data)
	else:
		processed_data, cluster_data = all_rxns(data)
	NMDS_data = perform_NMDS(processed_data, cluster_data)
	plot_pathogen(NMDS_data, cluster_data, bacteria)












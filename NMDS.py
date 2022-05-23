#import cobra.test
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
from sklearn.cluster import KMeans
import argparse
import natsort

#Parse the users input arguments
'''parser = argparse.ArgumentParser(description = 'Reads .sbml models, performs flux sampling, and plots flux data, clustering on a certian input (pathogen name, oxygen usage status, gram')
parser.add_argument('--cluster_type', default = 'none')
parser.add_argument('--sample_size', default = '500')
parser.add_argument('--bypass', default = 'false')
parser.add_argument('--consensus_rxns', default = 'false')
args = parser.parse_args()'''



def read_flux_data(downsample):
	isolation_source_data = pd.read_excel('isolationsourcechron.xlsx')
	class_data = pd.read_excel('classchron.xlsx')
	phylum_data = pd.read_excel('phylumchron.xlsx')
	reactions = pd.read_excel('reactions.xlsx')
	taxids = pd.read_excel('Taxids.xlsx')
	print

	columnlist = reactions.loc[:,'Reaction'].tolist().append(['Isolation Source', 'Class', 'Phylum','Taxid'])
	print(columnlist)
	filelist = glob.glob("*.csv")
	filelist = natsort.natsorted(filelist)#filelist.sort(key = lambda f: int(re.sub('\D', '', f)))
	

	NMDS_data = pd.DataFrame(columns = columnlist)


	i = 0
	for file in filelist:
		print(i)
		data = pd.read_csv(file)
		data = pd.DataFrame(data)
		data = data.iloc[:int(downsample)]
		data['Isolation Source'] = isolation_source_data.loc[i,'Isolation Source']
		data['Class'] = class_data.loc[i,'Class']
		data['Phylum'] = phylum_data.loc[i, 'Phylum']
		data['Taxid'] = taxids.loc[i,'Taxids']
		i+=1
		NMDS_data = pd.concat([NMDS_data,data], ignore_index = True).fillna(0)

	print(NMDS_data)
	NMDS_Data_isosource = NMDS_data[['Isolation Source']]
	del NMDS_data['Isolation Source']

	NMDS_Data_taxid = NMDS_data[['Taxid']]
	del NMDS_data['Taxid']

	NMDS_Data_class = NMDS_data[['Class']]
	del NMDS_data['Class']

	NMDS_Data_phylum = NMDS_data[['Phylum']]
	del NMDS_data['Phylum']

	print('performing nmds')
	nmds = MDS(n_components = 1, n_init = 1, metric = True, random_state = 0, dissimilarity = 'precomputed')
	similarities = pairwise_distances(NMDS_data, metric = 'braycurtis')
	data_trasnformed = nmds.fit_transform(similarities)

	df = pd.concat([NMDS_Data_class, NMDS_data], axis = 1)
	df = df.apply(lambda x: pd.Series(x.dropna().values)).fillna('')

	df.loc[:,0:1].to_csv('nmdstoplot.csv')

	return

def plot_nmds():
	df = pd.read_excel('nmdstoplot.xlsx')

	targets = ['Gammaproteobacteria', 'Betaproteobacteria', 'Epsilonproteobacteria', 'Alphaproteobacteria', 'Bacilli', 'Clostridia', 'Negativicutes', 'Tissierellia', 'Erysipelotrichia', 'Actinomycetia', 'Coriobacteria', 'Bacteroidia', 'Flavobacteriia', 'Spirochaetia', 'Mollicutes', 'Fusobacteriia', 'Chlamydia', '<not present>']
	legend_names =  ['Gammaproteobacteria', 'Betaproteobacteria', 'Epsilonproteobacteria', 'Alphaproteobacteria', 'Bacilli', 'Clostridia', 'Negativicutes', 'Tissierellia', 'Erysipelotrichia', 'Actinomycetia', 'Coriobacteria', 'Bacteroidia', 'Flavobacteriia', 'Spirochaetia', 'Mollicutes', 'Fusobacteriia', 'Chlamydia', '<not present>']
	colors = ['#609abf', '#6bc7db', '#67b5a4', '#4acfb2', '#60d168', '#a1eb86', '#8eb555', '#d2de5b', '#edd253', '#ffb84f', '#ff8c12', '#ff6e4f', '#d9455a', '#eb50a0', '#b15adb', '#6e47b5']

	plt.figure(figsize = (15,10))

	for target, color in zip(targets, colors):
		indicesToKeep = df['Class'] == target
		plt.scatter(df.loc[indicesToKeep,0], df.loc[indicesToKeep,1], c = color, alpha = 0.5)

	NMDS1 = df.loc[indicesToKeep,0]
	NMDS2 = df.loc[indicesToKeep,1]

	plt.legend(legend_names)

	plt.xlabel('NMDS1')
	plt.ylabel('NMDS2')

	plt.savefig('class.png', dpi = 300)

	return

if __name__ =='__main__':
	read_flux_data(10)
	plot_nmds()
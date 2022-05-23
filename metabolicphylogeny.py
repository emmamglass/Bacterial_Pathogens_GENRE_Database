import numpy as np 
import pandas as pd
import argparse
import cobra
import glob
import csv
from sklearn.neighbors import DistanceMetric
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from sklearn.metrics.pairwise import pairwise_distances
from kneed import KneeLocator
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

#parse the input arguments
parser = argparse.ArgumentParser(description = 'Reads .sbml models, determines reactions, calculates jaccard distance, performs nmds on reactions content')
parser.add_argument('--rxn_list', default = 'yes')
parser.add_argument('--rxn_presence', default = 'yes')
parser.add_argument('--jaccard_distance', default = 'yes')
parser.add_argument('--jaccard_heatmap', default = 'yes')
parser.add_argument('--nmds', default = 'yes')
parser.add_argument('--tree', default = 'yes' )
args = parser.parse_args()

#first we will read in all sbmls, grab rxns present
def make_rxn_list():
	file_list = glob.glob("*.sbml")
	rxn_list = []
	strain = 0
	for file in file_list:
		model = cobra.io.read_sbml_model(file, low_memory=False)
		rxns = model.reactions
		for rxn in rxns:
			reaction = str(rxn).split()
			reaction = str(reaction[0:1])
			reaction = str(reaction).replace("['", '')
			reaction = str(reaction).replace(":']",'')
			rxn_list.append(str(reaction))
		strain += 1
		print(strain)

	all_rxn_list = [i for n, i in enumerate(rxn_list) if i not in rxn_list[:n]]

	print(all_rxn_list)
	return all_rxn_list


def make_rxn_df():
	file_list = glob.glob("*.sbml")
	rxn_list= []
	strain = 0

	file = open('all_rxn_list.txt', 'r')
	reactome = file.read().split(',')
	reactome[-1] = reactome[-1].strip()

	rows = file_list
	dfrows = []
	for row in rows:
		row = row.strip('.sbml')
		dfrows.append(row)

	df = pd.DataFrame(index = dfrows, columns = reactome)

	row_num = 0
	for file in file_list:
		model = cobra.io.read_sbml_model(file, low_memory = False)
		reactions = model.reactions
		rxn_presence = []
		for i in range(len(reactome)):
			rxn = df.columns[i]
			if rxn in reactions:
				rxn_presence.append(1)
			else:
				rxn_presence.append(0)
		#print(rxn_presence)
		df.loc[dfrows[row_num]] = rxn_presence
		row_num+=1
		print(row_num)

	df.to_csv('reactionpresence.csv')
	return df

def jaccard(list1, list2):
	intersection = len(list(set(list1).intersection(list2)))
	union = (len(list1)+len(list2))-intersection
	return float(intersection/union)

def jaccard_distance():
	file_list = glob.glob("*.sbml")
	rows = file_list
	dfrows = []
	for row in rows:
		row = row.strip('.sbml')
		dfrows.append(row)

	df = pd.read_csv('reactionpresence.csv', low_memory = False)

	jaccdf = pd.DataFrame(index = dfrows, columns = dfrows )

	print(jaccdf)

	dist =  DistanceMetric.get_metric('hamming')

	jaccard_row = []
	row_num = 0
	for i in range(len(dfrows)):
		for j in range(len(dfrows)):
			list1 = pd.Series(df.loc[i])
			list2 = pd.Series(df.loc[j])
			jaccard_dist = abs(1-dist.pairwise([list1, list2]))
			jaccard_dist = jaccard_dist[0,1]
			jaccard_row.append(jaccard_dist)
		jaccdf.loc[dfrows[row_num]] = jaccard_row
		print(row_num)
		row_num+=1
		jaccard_row = []
	jaccdf.to_csv('jaccard_distances.csv')

def heatmap():
	dfjacc = pd.read_csv('jaccard_distances.csv', low_memory = False)

	#scaler = MinMaxScaler()
	#dfjacc = scaler.fit_transform(dfjacc)

	dfrxnpres = pd.read_csv('reactionpresence.csv', low_memory = False)

	#sns.heatmap(dfjacc, cmap = 'viridis', vmin = 0.7, vmax = 1, square = True, xticklabels = False, yticklabels = False)

	sns.clustermap(dfrxnpres, metric = 'hamming', cmap = 'copper', vmin = 0, xticklabels = False, yticklabels = False, vmax = 1, figsize = (15,8))

	plt.show()

def rxn_presence_nmds():
	df = pd.read_csv('reactionpresence.csv', low_memory = False)

	with open('colors.csv', newline='') as file:
		colors = file.read().splitlines()

	nmds = MDS(n_components = 2, n_init =1, metric = True, random_state = 0, dissimilarity = 'precomputed')

	similarities = pairwise_distances(df, metric = 'hamming')
	data_transformed = nmds.fit_transform(similarities)

	df = pd.DataFrame(data_transformed)
	
	df = df.apply(lambda x: pd.Series(x.dropna().values)).fillna('')

	targets = df.index
	#legend_names = ['E. c. 1', 'E. c. 2', 'E. c. 3', 'E. c. 4', 'E. c. 5', 'M. t. 1', 'M. t. 2', 'M. t. 3', 'M. t. 3', 'M. t. 4', 'M. t. 5', 'R. r. 1', 'R. r. 2', 'R. r. 3', 'R. r. 4', 'R. r. 5', 'P. a. 1', 'P. a. 2', 'P. a. 3', 'P. a. 4', 'P. a. 5', 'iPau21 a', 'iPau21 b', 'iPau21 c', 'iPau21 d', 'E. coli A', 'E. coli B', 'iPAU1129']
	
	colors = colors

	plt.figure(figsize=(15,8))
	
	#for target, color in zip(targets, colors):
	for target in targets:
		plt.scatter(df.loc[target,0],df.loc[target,1], c = '#1c5e9c', alpha = 0.7)

	NMDS1 = df.loc[target,0]
	NMDS2 = df.loc[target,1]

	#plt.legend(legend_names)

	plt.xlabel('NMDS1')
	plt.ylabel('NMDS2')

	#plt.title('NMDS')

	#plt.savefig('pathogen.png', dpi = 300)
	#plt.show()
	return data_transformed

def kmeans_clustering(clusters, data_transformed):
	kmeans = KMeans(n_clusters = int(clusters))

	kmeans.fit(data_transformed)
	y_kmeans = kmeans.predict(data_transformed)

	plt.scatter(data_transformed[:,0], data_transformed[:,1], c = y_kmeans, cmap = 'turbo', alpha = 0.7)
	centers = kmeans.cluster_centers_
	plt.scatter(centers[:,0], centers[:,1], c = '#000000', alpha = 0.1, s=1000)

	df = pd.DataFrame(y_kmeans)

	df.to_csv('y_kmeans.csv')
	plt.show()





#####################################################################################
# Main

#process input settings
get_rxn_list = str(args.rxn_list)
get_jaccard_distance = str(args.jaccard_distance)
rxn_presence = str(args.rxn_presence)
make_jaccard_heatmap = str(args.jaccard_heatmap)
do_nmds = str(args.nmds)
make_tree = str(args.tree)

if get_rxn_list == 'yes':
	rxns = make_rxn_list()
	with open('all_rxn_list', 'w') as f:
		write = csv.writer(f)
		write.writerow(rxns)

if rxn_presence == 'yes':
	make_rxn_df()

if get_jaccard_distance == 'yes':
	jaccard_distance()

if make_jaccard_heatmap == 'yes':
	heatmap()

if do_nmds == 'yes':
	data_transformed = rxn_presence_nmds()
	kmeans_clustering(3, data_transformed)






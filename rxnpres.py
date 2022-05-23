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
import matplotlib as mpl
import scipy.cluster.hierarchy as shc
import plotly.figure_factory as ff
from sklearn.manifold import TSNE

#parse the input arguments
parser = argparse.ArgumentParser(description = 'Reads .sbml models, determines reactions, calculates jaccard distance, performs nmds on reactions content')
parser.add_argument('--rxn_list', default = 'no')
parser.add_argument('--rxn_presence', default = 'no')
parser.add_argument('--jaccard_distance', default = 'no')
parser.add_argument('--jaccard_heatmap', default = 'no')
parser.add_argument('--nmds', default = 'no')
parser.add_argument('--tree', default = 'no' )
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

def make_gene_list():
	file_list = glob.glob("*.sbml")
	gene_list = []
	strain = 0
	for file in file_list:
		model = cobra.io.read_sbml_model(file, low_memory=False)
		genes = model.genes
		for gene in genes:
			gen = str(gene).split()
			gen = str(gen[0:1])
			gen = str(gen).replace("['", '')
			gen = str(gen).replace("']",'')
			gene_list.append(str(gen))
		strain += 1
		print(strain)

	all_gene_list = [i for n, i in enumerate(gene_list) if i not in gene_list[:n]]

	print(all_gene_list)
	return all_gene_list


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

def is_unique(s):
	a = s.to_numpy()
	return (a[0] == a).all()

def make_periph_metabol():

	df = pd.read_csv('reactionpresence.csv')
	print(df.shape)

	for column in df:
		if is_unique(df[column]) == True:
			df = df.drop(column, axis = 1)

	print(df.shape)


	df.to_csv('peripheral_rxns.csv')

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

	df = pd.read_csv('peripheral_rxns.csv', low_memory = False)

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
	jaccdf.to_csv('jaccard_distances_peripheral_rxns.csv')

def heatmap():
	dfjacc = pd.read_csv('jaccard_distances.csv', low_memory = False)

	#scaler = MinMaxScaler()
	#dfjacc = scaler.fit_transform(dfjacc)

	dfrxnpres = pd.read_csv('reactionpresence.csv', low_memory = False)

	cmap = mpl.colors.ListedColormap(['#112151', '#5FB861'])
	bounds = [0., 0.5, 1.]
	norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

	#sns.heatmap(dfjacc, cmap = 'viridis', vmin = 0.6, vmax = 1, square = True, xticklabels = False, yticklabels = False)

	#sns.clustermap(dfrxnpres, metric = 'hamming', cmap = cmap, vmin = 0, xticklabels = False, yticklabels = False, vmax = 1, figsize = (15,8), norm=norm)
    
	shc.dendrogram(shc.linkage(dfrxnpres, method='ward'))

	#fig = ff.create_dendrogram(dfrxnpres, orientation = 'left')

	#fig.show()

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


#This plots the results of the NMDS and clusters based upon pathogen name
def plot_isolation_source():
	isolation_source_data = pd.read_csv('isolationsourcedata.csv', low_memory = False)
	isolation_source_data = pd.DataFrame(isolation_source_data)

	df = pd.read_csv('reactionpresence1.csv', low_memory = False)
	df = pd.concat([isolation_source_data, df], axis = 1)
	df = df.apply(lambda x: pd.Series(x.dropna().values)).fillna('')
	#df = df.mask(df.eq('None')).dropna()
	#df = df.mask(df.eq('Stool')).dropna()
	#df = df.mask(df.eq('Blood')).dropna()
	'''df = df.mask(df.eq('Mouth')).dropna()
	df = df.mask(df.eq('Lymph Node')).dropna()
	df = df.mask(df.eq('Urine')).dropna()
	df = df.mask(df.eq('Pus')).dropna()
	df = df.mask(df.eq('Abdomen')).dropna()
	df = df.mask(df.eq('Vagina')).dropna()
	df = df.mask(df.eq('Liver')).dropna()
	df = df.mask(df.eq('Semen')).dropna()
	df = df.mask(df.eq('Kidney')).dropna()
	df = df.mask(df.eq('Ear')).dropna()
	df = df.mask(df.eq('Eye')).dropna()
	df = df.mask(df.eq('Axillia')).dropna()
	df = df.mask(df.eq('Bone')).dropna()
	df = df.mask(df.eq('Joint')).dropna()
	df = df.mask(df.eq('Skin')).dropna()
	df = df.mask(df.eq('Abscess')).dropna()
	df = df.mask(df.eq('Wound')).dropna()
	df = df.mask(df.eq('Burn')).dropna()
	df = df.mask(df.eq('Brain')).dropna()
	df = df.mask(df.eq('Cerebrospinal Fluid')).dropna()'''
	new_iso = df['Isolation Source'].values
	new_iso = pd.DataFrame(new_iso)
	new_iso.columns = ['Isolation Source']
	print(new_iso)
	df = df.drop('Isolation Source', axis = 1)

	#nmds = MDS(n_components = 2, n_init =1, metric = True, random_state = 0, dissimilarity = 'precomputed')
	nmds = TSNE(n_components = 2, perplexity=30.0, early_exaggeration=12.0, learning_rate=200.0,  n_iter=1000, init='random', method='exact')

	similarities = pairwise_distances(df, metric = 'hamming')

	data_transformed = nmds.fit_transform(similarities)

	data = data_transformed

	df = pd.DataFrame(data)
	df = pd.concat([new_iso, df], axis = 1)
	print(df)

	targets = ['None',
			   'Lung', 'Respiratory', 'Sputum','Throat',
			   'Gut', 'Stomach', 'Bile', 
			   'Brain', 'Cerebrospinal Fluid', 
			   'Bone', 'Joint',
			   'Skin', 'Abscess','Wound', 'Burn',
			   'Blood', 'Stool', 'Mouth', 'Lymph Node', 'Urine', 'Pus', 'Abdomen', 'Vagina', 'Liver','Semen', 'Kidney', 'Ear', 'Eye', 'Axillia']
	#targets = ['Non-CF', 'CF']
	#legend_names = ['Non-CF', 'CF']
	#colors = ['#18851a', '#14a395']
	legend_names =  ['None',
					 'Lung', 'Respiratory', 'Sputum','Throat',
			   		 'Gut', 'Stomach', 'Bile', 
			   		 'Brain', 'Cerebrospinal Fluid', 
			   		 'Bone', 'Joint',
			   		 'Skin', 'Abscess','Wound', 'Burn',
			   		 'Blood', 'Stool', 'Mouth', 'Lymph Node', 'Urine', 'Pus', 'Abdomen', 'Vagina', 'Liver','Semen', 'Kidney', 'Ear', 'Eye', 'Axillia']
	colors = ['#dedede',
			  '#0388fc', '#005dad', '#014885', '#00335e', 
			  '#8fcf5b', '#6ea343', '#46692b',
			  '#8153b5', '#4f2d75',
			  '#edce00', '#baa200',
			  '#c94e24', '#ad411c', '#8f3111', '#782205',
			  '#bd0000', '#4A1C00', '#edd18e', '#ca8fcc', '#e8e351', '#5ebf80', '#5ebfb7', '#ff9063', '#513f8a', '#949494', '#95aba2', '#747846', '#856e8a', '#fceb95']
	plt.figure(figsize=(15,8))
	
	for target, color in zip(targets, colors):
		indicesToKeep = df['Isolation Source'] == target
		#print(indicesToKeep)
		plt.scatter(df.loc[indicesToKeep,1],-df.loc[indicesToKeep,0], c = color, alpha = .9, s = 50)

	NMDS1 = df.loc[indicesToKeep,0]
	NMDS2 = df.loc[indicesToKeep,1]

	#plt.legend(legend_names, fontsize = 8)

	#plt.xlabel('NMDS1')
	#plt.ylabel('NMDS2')

	plt.xlabel('tSNE1')
	plt.ylabel('tSNE2')


	#plt.title('NMDS')

	#plt.savefig('pathogen.png', dpi = 300)
	plt.show()
	return

def plot_phylum():
	phylum_data = pd.read_csv('phylum.csv', low_memory = False)
	phylum_data = pd.DataFrame(phylum_data)

	df = pd.read_csv('reactionpresence1.csv', low_memory = False)
	df = pd.concat([phylum_data, df], axis = 1)
	df = df.apply(lambda x: pd.Series(x.dropna().values)).fillna('')

	new_phy = df['Phylum'].values
	new_phy = pd.DataFrame(new_phy)
	new_phy.columns = ['Phylum']
	print(new_phy)
	df = df.drop('Phylum', axis = 1)


	nmds = TSNE(n_components = 2, perplexity=30.0, early_exaggeration=12.0, learning_rate=200.0,  n_iter=1000, init='random', method='exact')

	#nmds = MDS(n_components = 2, n_init =1, metric = True, random_state = 0, dissimilarity = 'precomputed')

	similarities = pairwise_distances(df, metric = 'hamming')

	data_transformed = nmds.fit_transform(similarities)

	data = data_transformed

	df = pd.DataFrame(data)
	df = pd.concat([new_phy, df], axis = 1)
	print(df)

	targets = ['Actinobacteria', 'Bacteroidetes', 'Chlamydiae', 'Firmicutes', 'Fusobacteria', 'Proteobacteria', 'Saccharibacteria', 'Spirochaetes', 'Tenericutes']

	#targets = ['Non-CF', 'CF']
	#legend_names = ['Non-CF', 'CF']
	#colors = ['#18851a', '#14a395']
	legend_names =  ['Actinobacteria', 'Bacteroidetes', 'Chlamydiae', 'Firmicutes', 'Fusobacteria', 'Proteobacteria', 'Saccharibacteria', 'Spirochaetes', 'Tenericutes']
	
	colors = ['#ca67a0', '#d0413d', '#ef8376', '#ed7720', '#ffbc41', '#83996b', '#619299', '#004c63', '#0e204e']
	
	plt.figure(figsize=(15,8))
	
	for target, color in zip(targets, colors):
		indicesToKeep = df['Phylum'] == target
		#print(indicesToKeep)
		plt.scatter(df.loc[indicesToKeep,1],-df.loc[indicesToKeep,0], c = color, alpha = 0.9, s = 50)

	NMDS1 = df.loc[indicesToKeep,0]
	NMDS2 = df.loc[indicesToKeep,1]

	plt.legend(legend_names, fontsize = 8)

	#plt.xlabel('NMDS1')
	#plt.ylabel('NMDS2')

	plt.xlabel('tSNE1')
	plt.ylabel('tSNE2')

	#plt.title('NMDS')

	#plt.savefig('pathogen.png', dpi = 300)
	plt.show()
	return

def plot_class():
	class_data = pd.read_csv('class.csv', low_memory = False)
	class_data = pd.DataFrame(class_data)

	df = pd.read_csv('reactionpresence1.csv', low_memory = False)
	df = pd.concat([class_data, df], axis = 1)
	df = df.apply(lambda x: pd.Series(x.dropna().values)).fillna('')

	new_class = df['Class'].values
	new_class = pd.DataFrame(new_class)
	new_class.columns = ['Class']
	print(new_class)
	df = df.drop('Class', axis = 1)

	#nmds = MDS(n_components = 2, n_init =1, metric = True, random_state = 0, dissimilarity = 'precomputed')

	nmds = TSNE(n_components = 2, perplexity=30.0, early_exaggeration=12.0, learning_rate=200.0,  n_iter=1000, init='random', method='exact')

	similarities = pairwise_distances(df, metric = 'hamming')

	data_transformed = nmds.fit_transform(similarities)

	data = data_transformed

	df = pd.DataFrame(data)
	df = pd.concat([new_class, df], axis = 1)
	print(df)

	targets = ['Gammaproteobacteria', 'Betaproteobacteria', 'Epsilonproteobacteria', 'Alphaproteobacteria', 'Bacilli', 'Clostridia', 'Negativicutes', 'Tissierellia', 'Erysipelotrichia', 'Actinomycetia', 'Coriobacteria', 'Bacteroidia', 'Flavobacteriia', 'Spirochaetia', 'Mollicutes', 'Fusobacteriia', 'Chlamydia', '<not present>']

	#targets = ['Non-CF', 'CF']
	#legend_names = ['Non-CF', 'CF']
	#colors = ['#18851a', '#14a395']
	legend_names =  ['Gammaproteobacteria', 'Betaproteobacteria', 'Epsilonproteobacteria', 'Alphaproteobacteria', 'Bacilli', 'Clostridia', 'Negativicutes', 'Tissierellia', 'Erysipelotrichia', 'Actinomycetia', 'Coriobacteria', 'Bacteroidia', 'Flavobacteriia', 'Spirochaetia', 'Mollicutes', 'Fusobacteriia', 'Chlamydia', '<not present>']
	
	colors = ['#609abf', '#6bc7db', '#67b5a4', '#4acfb2', '#60d168', '#a1eb86', '#8eb555', '#d2de5b', '#edd253', '#ffb84f', '#ff8c12', '#ff6e4f', '#d9455a', '#eb50a0', '#b15adb', '#6e47b5']

	plt.figure(figsize=(15,8))
	
	for target, color in zip(targets, colors):
		indicesToKeep = df['Class'] == target
		#print(indicesToKeep)
		plt.scatter(df.loc[indicesToKeep,1],-df.loc[indicesToKeep,0], c = color, alpha = .9, s = 50)

	NMDS1 = df.loc[indicesToKeep,0]
	NMDS2 = df.loc[indicesToKeep,1]

	plt.legend(legend_names, fontsize = 8)

	plt.xlabel('NMDS1')
	plt.ylabel('NMDS2')

	plt.xlabel('tSNE1')
	plt.ylabel('tSNE2')

	#plt.title('NMDS')

	#plt.savefig('pathogen.png', dpi = 300)
	plt.show()
	return



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

#make_periph_metabol()
#jaccard_distance()
#plot_isolation_source()
#plot_phylum()
#plot_class()
plot_class()







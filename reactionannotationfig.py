import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

#load reaction presence data
rxnpres_data = pd.read_csv('reactionpresence.csv')
#remove the _c from reactions 
rxnpres_data.columns = rxnpres_data.columns.str.replace('_c', '')
#sum
rxnpres_data.loc['total'] = rxnpres_data.sum()
#add reaction counts
counts = rxnpres_data.loc['total']
counts = counts[1:]
#grab reactions
countsrxns = rxnpres_data.columns
countsrxns = countsrxns[1:]
#make counts/reactions df
Countsdf = pd.DataFrame([countsrxns, counts]).T
Countsdf.columns = ['Reaction', 'Counts']

#load annotation data
reaction_annotation = pd.read_csv('reaction_annotations.csv')
#grab annotation column and reaction column
annotation = reaction_annotation.loc[:,'Annotation']
annotationrxns = reaction_annotation.loc[:, 'Unnamed:']
#make annotation/reactions dataframe
Annotationdf = pd.DataFrame([annotationrxns, annotation]).T
Annotationdf.columns = ['Reaction', 'Annotation']

#merge
annotatedrxns = pd.merge(Countsdf, Annotationdf, on='Reaction', how='left')
annotatedrxns = annotatedrxns.replace({np.nan: 'None'})
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.strip("'>")
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.strip('n')
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.strip('\\')
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.replace("match='Metabolism; Carbohydrate metabolism", 'Carbohydrate Metabolism')
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.replace("match='Metabolism; Energy metabolism", 'Energy Metabolism')
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.replace("match='Metabolism; Amino acid metabolism", 'Amino Acid Metabolism')
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.replace("match='Metabolism; Xenobiotics biodegradation and metabo", 'Xenobiotics Degredation and Metabolism')
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.replace("match='Metabolism; Lipid metabolism", 'Lipid Metabolism')
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.replace("match='Metabolism; Metabolism of cofactors and vitamins", 'Metabolism of Cofactors and Vitamins')
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.replace("match='Metabolism; Metabolism of other amino acids", 'Metabolism of Other Amino Acids')
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.replace("match='Metabolism; Glycan biosynthesis and metabolism", 'Glycan Biosynthesis and Metabolism')
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.replace("match='Metabolism; Biosynthesis of other secondary metab", 'Biosynthesis of Other Secondary Metabolites')
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.replace("match='Metabolism; Metabolism of terpenoids and polyketi", 'Metabolism of Terpenoids and Polyketides')
annotatedrxns['Annotation'] = annotatedrxns['Annotation'].str.replace("match='Metabolism; Nucleotide metabolism", 'Nucleotide Metabolism')
annotatedrxns = annotatedrxns[annotatedrxns['Annotation'] != 'None']


'''Carbohydrate = annotatedrxnstoplot.loc[annotatedrxnstoplot['Annotation'] == "Carbohydrate Metabolism"]
Energy = annotatedrxns.loc[annotatedrxns['Annotation'] == "Energy Metabolism"]
#print(Energy)
nan = annotatedrxns.loc[annotatedrxns['Annotation'] == "None"]
#print(nan)
AminoAcid = annotatedrxns.loc[annotatedrxns['Annotation'] == "match='Metabolism; Amino acid metabolism'>"]
Nucleotide = annotatedrxns.loc[annotatedrxns['Annotation'] == "match='Metabolism; Nucleotide metabolism'>"]
Xenobiotics = annotatedrxns.loc[annotatedrxns['Annotation'] == "match='Metabolism; Xenobiotics biodegradation and metabo'>"]
Lipid = annotatedrxns.loc[annotatedrxns['Annotation'] == "match='Metabolism; Lipid metabolism'>"]
Cofactors = annotatedrxns.loc[annotatedrxns['Annotation'] == "match='Metabolism; Metabolism of cofactors and vitamins'>"]
OtherAminoAcid = annotatedrxns.loc[annotatedrxns['Annotation'] == "match='Metabolism; Metabolism of other amino acids'>"]
Glycan = annotatedrxns.loc[annotatedrxns['Annotation'] == "match='Metabolism; Glycan biosynthesis and metabolism'>"]
Secondary =annotatedrxns.loc[annotatedrxns['Annotation'] == "match='Metabolism; Biosynthesis of other secondary metab'>"]
Terpenoids = annotatedrxns.loc[annotatedrxns['Annotation'] == "match='Metabolism; Metabolism of terpenoids and polyketi'>"]

rxnpres_data.loc['annotation'] = annotation'''


n_bins = 8

x = counts[1:]

fig, axes= plt.subplots(figsize = (9,8))

#colors = [ '#2a307a', '#6a5b9c', '#9c9bd5', '#b07798', '#f7b4b1', '#e5806e', '#c03323', '#e69a00', '#f6bc4e', '#ada236', '#315724']#'#13124a',
#colors = ['#801012', '#c71e1b', '#f15a3f', '#fa8261', '#ffad7f', '#f9dfa4', '#ced2b1', '#97ab75', '#738d43', '#436027', '#224b20']
colors = ['#5f9bbf', '#68b26b', '#7db66a', '#9bba6f', '#bbc462', '#bbc462', '#fcad59', '#e77961', '#d86878', '#a776bc', '#9780c0' ]

x = annotatedrxns.pivot(columns='Annotation')['Counts']

axes.hist(x, n_bins, histtype='bar', stacked=True, color = colors, edgecolor = 'black')

axes.legend(['Amino Acid Metabolism', 'Biosynthesis of Other Secondary Metabolites', 'Carbohydrate Metabolism', 'Energy Metabolism', 'Glycan Biosynthesis and Metabolism',
			'Lipid Metabolism', 'Metabolism of Cofactors and Vitamins', 'Metabolism of Other Amino Acids', 'Metabolism of Terpenoids and Polyketides', 
			'Nucleotide Metabolism', 'Xenobiotics Degredation and Metabolism'])#, 'None'])

#axes.set_yscale('log')
axes.set_xlabel('Number of Models Containing a Reaction')
axes.set_ylabel('Number of Reactions')



plt.show()

#func = pd.concat([counts,annotation], axis = 1, ignore_index= True)
#print(func)

'''uniquefunction = ["match='Metabolism; Carbohydrate metabolism\n'>", "match='Metabolism; Energy metabolism\n'>", 
				  "", "match='Metabolism; Amino acid metabolism\n'>", "match='Metabolism; Nucleotide metabolism\n'>",
				  "match='Metabolism; Xenobiotics biodegradation and metabo\n'>", "match='Metabolism; Lipid metabolism\n'>",
				  "match='Metabolism; Metabolism of cofactors and vitamins\n'>", "match='Metabolism; Metabolism of other amino acids\n'>",
				  "match='Metabolism; Glycan biosynthesis and metabolism\n'>", "match='Metabolism; Biosynthesis of other secondary metab\n'>",
				  "match='Metabolism; Metabolism of terpenoids and polyketi\n'>"]

i = 0
histogramdata = pd.DataFrame(columns = ['Carbohydrate Metabolism', 'Energy Metabolism', '', 'Amino Acid Metabolism',
	    								'Nucleotide Metabolism', 'Xenobiotics Biodegradation and Metabolism',
	    								'Lipid Metabolism', 'Metabolism of Cofactors and Vitamins', 'Metabolism of Other Amino Acids', 
	    								'Glycan Biosynthesis and Metabolism', 'Biosynthesis of Other Secondary Metabolites',
	    								'Metabolism of Terpenoids and Polyketides'])



n_bins = 30

x = counts[1:]

fig, axes= plt.subplots(figsize = (9,8))

colors = ['#012A4A', '#013A63', '#01497C', '#014F86', '#2A6F97', '#2C7DA0', '#468FAF', '#61A5C2', '#89C2D9', '#A9D6E5']

axes.hist(x, n_bins, histtype='bar', stacked=True, color = '#2C7DA0', edgecolor = 'black')
axes.set_yscale('log')
axes.set_xlabel('Number of Models Containing a Reaction')
axes.set_ylabel('Number of Reactions')'''



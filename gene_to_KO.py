from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
import pandas as pd
import numpy as np
import math



gene_essentiality = pd.read_csv('EssentialGenes.csv')

df = pd.DataFrame()
for column in gene_essentiality.columns:
	gene_names = []
	for gene in gene_essentiality[column]:
		if pd.isnull(gene) == False:
			print(gene)
			try:
				kegg = REST.kegg_get(gene).read()
				position = kegg.find("ORTHOLOGY")
				if position != -1:
					gene_names.append(kegg[position+12:position+18])
				elif position == -1:
					gene_names.append("NA")
			except:
				gene_names.append('No KEGG Ortholog Found')
				print('WE HAVE ENTERED THE EXCEPT BLOCK')
	df[column] = pd.Series(gene_names)
print(df)

df.to_csv('KO_essentialgenes.csv')

#gene_essentiality.index = gene_names
#gene_essentiality.to_csv("svm_gene_essentiality.csv", header = list(essential_genes_df.columns))
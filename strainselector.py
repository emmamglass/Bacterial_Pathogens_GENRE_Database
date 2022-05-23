import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

#read in dataframe
df = pd.read_excel("All_Bacterial_Pathogen_Strains.xlsx", dtype = {'Genome ID': np.unicode_})
print(df)

#determine each unique taxon id in the larger dataset and create a vector
taxid = df['NCBI Taxon ID']
uniquetaxid = taxid.unique()

#how many unique taxon ids there are....
print('There are ' + str(len(uniquetaxid)) + ' unique NCBI Taxon IDs.')

#Lets write a program that chooses one strain per NCBI Taxon ID
#We will select strains that have the most metadata associated with them

master = pd.DataFrame()

num_strains_tot = []
num_index = []


#for each unique taxon ID, pull out all strains with corresponding taxon ID
progress = 0
for num in uniquetaxid:
	temp = pd.DataFrame()
	
	for index, row in df.iterrows():
		if str(num) == str(df.iloc[index]['NCBI Taxon ID']):
			temp = temp.append(row.transpose())
			df.drop(labels = index, axis = 0)
	progress += 1
	print(str(int((progress/914)*100)) + '%')

	count = temp.count(axis = 1)
	
	#Determine how many strains there are per ncbi taxon id and save that to a dictionary
	num_strains = temp.shape[0]
	num_strains_tot.append(int(num_strains))

	#Choosing strain with the most metadata, ranking based on metadata importance
	for index, row in temp.iterrows():
		new_count = 0
		new_row = pd.DataFrame()
		if (count[index] > new_count):
			new_count = count[index]
			new_row = row 
			if ((pd.isnull(df.iloc[index]['Isolation Source'])) == False):
				new_count = count[index]
				new_row = row 
				if ((pd.isnull(df.iloc[index]['Host Health'])) == False):
					new_count = count[index]
					new_row = row 
					if ((pd.isnull(df.iloc[index]['Isolation Country'])) ==False):
						new_count = count[index]
						new_row = row 
						if ((pd.isnull(df.iloc[index]['Collection Date']))==False):
							new_count = count[index]
							new_row = row
							if ((pd.isnull(df.iloc[index]['Host Age'])) ==False):
								new_count = count[index]
								new_row = row 
	#add chosen strain to masterdoc				
	master = master.append(new_row)

master = master.convert_dtypes(convert_string = True)
print(master)
counts = []


#turn masterdoc into a sv
master.to_csv('Selected_Strains.csv', float_format = '{:}')



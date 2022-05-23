
import pandas as pd

with open('Group1_k20.csv', newline = '') as file:
	Group1 = file.read().splitlines()
with open('Group2_k20.csv', newline = '') as file:
	Group2 = file.read().splitlines()
with open('Group3_k20.csv', newline = '') as file:
	Group3 = file.read().splitlines()
with open('Group4_k20.csv', newline = '') as file:
	Group4 = file.read().splitlines()
with open('Group5_k20.csv', newline = '') as file:
	Group5 = file.read().splitlines()
with open('Group6_k20.csv', newline = '') as file:
	Group6 = file.read().splitlines()
with open('Group7_k20.csv', newline = '') as file:
	Group7 = file.read().splitlines()
with open('Group8_k20.csv', newline = '') as file:
	Group8 = file.read().splitlines()
with open('Group9_k20.csv', newline = '') as file:
	Group9 = file.read().splitlines()
with open('Group10_k20.csv', newline = '') as file:
	Group10 = file.read().splitlines()
with open('Group11_k20.csv', newline = '') as file:
	Group11 = file.read().splitlines()
with open('Group12_k20.csv', newline = '') as file:
	Group12 = file.read().splitlines()
with open('Group13_k20.csv', newline = '') as file:
	Group13 = file.read().splitlines()
with open('Group14_k20.csv', newline = '') as file:
	Group14 = file.read().splitlines()
with open('Group15_k20.csv', newline = '') as file:
	Group15 = file.read().splitlines()
with open('Group16_k20.csv', newline = '') as file:
	Group16 = file.read().splitlines()
with open('Group17_k20.csv', newline = '') as file:
	Group17 = file.read().splitlines()
with open('Group18_k20.csv', newline = '') as file:
	Group18 = file.read().splitlines()
with open('Group19_k20.csv', newline = '') as file:
	Group19 = file.read().splitlines()

df = pd.read_csv('Selected_Strains.csv')

Group1arr = pd.DataFrame()

i = 0
for taxid in Group1:
	for row in df.iterrows():
		if taxid == df.iloc([i,'NCBI Taxon ID']):
			print(taxid)



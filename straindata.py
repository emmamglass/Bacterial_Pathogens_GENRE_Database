import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

#read in dataframe
df = pd.read_csv("Selected_Strains.csv")

isolation_country = df['Isolation Country'].tolist()
host_health = df['Host Health'].tolist()
isolation_source = df['Isolation Source'].tolist()
cell_shape = df['Cell Shape'].tolist()

unique_isolation_country = df['Isolation Country'].unique().tolist()
unique_host_health = df['Host Health'].unique().tolist()
unique_isolation_source = df['Isolation Source'].unique().tolist()
unique_cell_shape = df['Cell Shape'].unique().tolist()


#print(len(unique_isolation_source))
#print('---------------')
#print(len(unique_host_health))
#print('---------------')
#print(len(unique_isolation_country))

#### isolation Country ######
count_list = []
for country in unique_isolation_country:
	counts = isolation_country.count(country)
	count_list.append(counts)

country_counts = {k:v for k,v in zip(unique_isolation_country, count_list)}
#country_counts = sorted(country_counts.items(), key=lambda x: x[1], reverse=True)
l = len(country_counts)
marklist=list(country_counts.items())
for i in range(l-1):
	for j in range(i+1, l):
		if marklist[i][1]>marklist[j][1]:
			t=marklist[i]
			marklist[i]=marklist[j]
			marklist[j]=t 
	country_counts = dict(marklist)

country_keys = country_counts.keys()
country_values = country_counts.values()

##### isolation source #####
count_list = []
for source in unique_isolation_source:
	counts = isolation_source.count(source)
	count_list.append(counts)

source_counts = {k:v for k,v in zip(unique_isolation_source, count_list)}
#country_counts = sorted(country_counts.items(), key=lambda x: x[1], reverse=True)
l = len(source_counts)
marklist=list(source_counts.items())
for i in range(l-1):
	for j in range(i+1, l):
		if marklist[i][1]>marklist[j][1]:
			t=marklist[i]
			marklist[i]=marklist[j]
			marklist[j]=t 
	source_counts = dict(marklist)

source_keys = source_counts.keys()
source_values = source_counts.values()


####### Host health ######
count_list = []
for host in unique_host_health:
	counts = host_health.count(host)
	count_list.append(counts)

host_counts = {k:v for k,v in zip(unique_host_health, count_list)}
#country_counts = sorted(country_counts.items(), key=lambda x: x[1], reverse=True)
l = len(host_counts)
marklist=list(host_counts.items())
for i in range(l-1):
	for j in range(i+1, l):
		if marklist[i][1]>marklist[j][1]:
			t=marklist[i]
			marklist[i]=marklist[j]
			marklist[j]=t 
	host_counts = dict(marklist)

host_keys = host_counts.keys()
host_values = host_counts.values()


##### Cell Shape #####



count_list = []
for shape in unique_cell_shape:
	counts = cell_shape.count(shape)
	count_list.append(counts)

cell_counts = {k:v for k,v in zip(unique_cell_shape, count_list)}
#country_counts = sorted(country_counts.items(), key=lambda x: x[1], reverse=True)
l = len(cell_counts)
marklist=list(cell_counts.items())
for i in range(l-1):
	for j in range(i+1, l):
		if marklist[i][1]>marklist[j][1]:
			t=marklist[i]
			marklist[i]=marklist[j]
			marklist[j]=t 
	cell_counts = dict(marklist)

cell_keys = cell_counts.keys()
cell_values = cell_counts.values()

country_counts.popitem()
source_counts.popitem()
host_counts.popitem()
cell_counts.popitem()

print(len(country_counts))
print('___________________')
print(len(source_counts))
print('___________________')
print(len(host_counts))
print('___________________')
print(len(cell_counts))



fig, ax = plt.subplots(1,1, figsize = (15, 7))
plt.barh(range(len(country_counts)), list(country_values),align = 'center', color = 'darkgreen')
plt.yticks(range(len(country_counts)), list(country_keys), rotation = 0, fontsize = 8)
plt.xscale('log')
plt.show()

fig, ax = plt.subplots(1,1, figsize = (15, 7))
plt.barh(range(50), list(source_values)[-50:],align = 'center', color = 'mediumseagreen')
plt.yticks(range(50), list(source_keys)[-50:], rotation = 0, fontsize = 8)
plt.xscale('log')
plt.show()

fig, ax = plt.subplots(1,1, figsize = (15, 7))
plt.barh(range(50), list(host_values)[-50:],align = 'center', color = 'lightgreen')
plt.yticks(range(50), list(host_keys)[-50:], rotation = 0, fontsize = 8)
plt.xscale('log')
plt.show()

fig, ax = plt.subplots(1,1, figsize = (15, 7))
plt.barh(range(len(cell_counts)), list(cell_values),align = 'center', color = 'darkolivegreen')
plt.yticks(range(len(cell_counts)), list(cell_keys), rotation = 0, fontsize = 8)
plt.xscale('log')
plt.show()



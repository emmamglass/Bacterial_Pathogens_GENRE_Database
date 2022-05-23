import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

class_data = pd.read_excel('nmds_class_boxplot.xlsx')
isosource_data = pd.read_excel('nmds_isosource_boxplot.xlsx')

fig, ax = plt.subplots(figsize = (15,8))

ax = sns.boxplot(data = class_data)
ax = sns.stripplot(data = class_data, color = 'k', size = 2)
plt.xticks(rotation = 30, ha = 'right')
plt.xlabel('Class', fontsize = 13, weight = 'bold')
plt.ylabel('Flux', fontsize = 13, weight = 'bold')

fig, ax = plt.subplots(figsize = (15,8))
ax = sns.boxplot(data = isosource_data)
ax = sns.stripplot(data = isosource_data, color = 'k', size = 2)
plt.xticks(rotation = 30, ha='right')
plt.xlabel('Isolation Source', fontsize = 13, weight = 'bold')
plt.ylabel('Flux', fontsize = 13, weight = 'bold')

plt.show()
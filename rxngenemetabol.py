import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

data = pd.read_csv('taxonomyinfo.csv')

fig, axes = plt.subplots(3,1, figsize = (9,8))

axes[0].axhline(y = 1407, color = '#B8B8B8', linestyle = ':', zorder = 0)
sns.boxplot(ax = axes[0], x='Phylum', y='Reactions', data = data, color = '#9bc283')
axes[0].set_xticklabels(['','','','','','','','',''])
axes[0].set_xlabel('')
axes[0].tick_params(axis = 'x', which = 'both', bottom = False)


axes[1].axhline(y = 1365, color = '#B8B8B8', linestyle = ':', zorder = 0)
sns.boxplot(ax=axes[1], x='Phylum', y = 'Genes', data = data, color = '#649c5f')
axes[1].set_xticklabels(['','','','','','','','',''])
axes[1].set_xlabel('')
axes[1].tick_params(axis = 'x', which = 'both', bottom = False)


axes[2].axhline(y = 1470, color = '#B8B8B8', linestyle = ':', zorder = 0)
sns.boxplot(ax = axes[2], x='Phylum', y = 'Metabolites', data = data, color = '#1a5920')
labels = ['Spirochaetes (27)', 'Proteobacteria (457)', 'Chlamydiae (9)', 'Bacteroidetes (39)', 'Fusobacteria (19)', 'Firmicutes (175)', 'Actinobacteria (168)', 'Tenericutes (19)', 'Saccharibacteria (1)']
axes[2].set_xticklabels(labels,rotation = 25)
axes[2].set_xlabel('')


plt.show()
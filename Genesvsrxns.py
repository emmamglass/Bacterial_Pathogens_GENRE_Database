import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

data = pd.read_csv('taxonomyinfo.csv')

fig, axes = plt.subplots(figsize = (3,2))

#axes.plot(y = 1407, color = '#B8B8B8', linestyle = ':', zorder = 0)
sns.scatterplot(ax = axes, x='Genes', y='Reactions', data = data, color = '#649c5f')
#axes.set_xticklabels(['','','','','','','','',''])
axes.set_xlabel('Number of Genes')
axes.set_ylabel('Number of Reactions')
#axes.tick_params(axis = 'x', which = 'both', bottom = False)

plt.show()
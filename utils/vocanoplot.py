%matplotlib inline
import pandas as pd # Data analysis
import numpy as np # Scientific computing
import matplotlib.pyplot as plt # Plotting
import matplotlib.colors as colors # Coloring
import seaborn as sns # Statistical visualization
fold1 = pd.DataFrame(fold)
fold1['a'] = range(len(fold1))
fold1 = fold1.set_index('a')
# In[*]
pvalue1 = pd.DataFrame(pvalue)

result = pd.concat([fold1, pvalue1],axis=1,ignore_index=True)  

result.columns = ['fold','pvalue']

result['log(pvalue)'] = -np.log10(result['pvalue'])

# In[*]
result['sig'] = 'normal'

result['size']  =np.abs(result['fold'])/10
 
result.loc[(result.fold> 1 )&(result.pvalue < 0.05),'sig'] = 'up'
result.loc[(result.fold< -1 )&(result.pvalue < 0.05),'sig'] = 'down'

# In[*]
 ax = sns.scatterplot(x="fold", y="log(pvalue)",
                      hue='sig',
                      hue_order = ('down','normal','up'),
                      palette=("#377EB8","grey","#E41A1C"),
                      data=result)
ax.set_ylabel('-log(pvalue)',fontweight='bold')
ax.set_xlabel('FoldChange',fontweight='bold')


fig, ax = plt.subplots(figsize = (6,6))


sc = sns.scatterplot(data = df.dropna(subset=['padj']), 
                    x = 'log2FoldChange', 
                    y = 'nlog10',
                    hue = 'updown', 
                    hue_order = ['NS','UP Genes', 'Down Genes'],
                    palette = ['lightgrey', 'orange', 'purple'],
                    #style = 'shape', 
                    #style_order = ['picked3', 'picked4', 'not_important'],
                    #markers = ['^', 's', 'o'], 
                    ax = ax,     
                    size = 'baseMean', sizes = (40, 400))

ax.axhline(2, zorder = 0, c = 'k', lw = 2, ls = '--')
ax.axvline(1, zorder = 0, c = 'k', lw = 2, ls = '--')
ax.axvline(-1, zorder = 0, c = 'k', lw = 2, ls = '--')


## add text
texts = []
for i in range(len(df)):
    if df.iloc[i].nlog10 > 20 and abs(df.iloc[i].log2FoldChange) > 4:
        texts.append(plt.text(x = df.iloc[i].log2FoldChange, y = df.iloc[i].nlog10, s = df.iloc[i].gene_name,
                             fontsize = 12, weight = 'bold'))
        
adjust_text(texts, arrowprops = dict(arrowstyle = '-', color = 'k'))





plt.legend(loc = 1, bbox_to_anchor = (1.4,1), frameon = False, prop = {'weight':'bold'})

# for axis in ['bottom', 'left']:
#     ax.spines[axis].set_linewidth(2)
    
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)

# ax.tick_params(width = 2)

# plt.xticks(size = 12, weight = 'bold')
# plt.yticks(size = 12, weight = 'bold')

# plt.xlabel("$log_{2}$ fold change", size = 15)
# plt.ylabel("-$log_{10}$ FDR", size = 15)

#plt.savefig('volcano.png', dpi = 300, bbox_inches = 'tight', facecolor = 'white')

plt.show()
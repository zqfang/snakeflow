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

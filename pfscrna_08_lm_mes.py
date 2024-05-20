import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager as fm, rcParams
##########
mes_pseudo = pd.read_csv('pseudo_mes.csv',header = 0, index_col=0)
sns.set_style('white')
sns.lmplot(scatter=False,data=mes_pseudo,x='Pseudotime',y='Gene_expression',
           hue='genes',
           #palette = neu_subset_palette,
           order=2,        
          )
plt.savefig('lm_mes_gene.pdf')
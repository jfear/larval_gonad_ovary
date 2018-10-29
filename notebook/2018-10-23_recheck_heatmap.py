import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

zscores = pd.read_parquet('../output/scrnaseq-wf/zscore_tpm.res.0.4.parquet')

#%%
sns.clustermap(zscores, cmap='viridis')
plt.show()

#%%
plt.plot(np.random.random(10))
plt.show()

#%%
plt.plot(np.random.random(10))
plt.show()
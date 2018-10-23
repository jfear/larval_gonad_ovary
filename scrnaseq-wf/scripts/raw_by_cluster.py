import pandas as pd

raw = pd.read_parquet(snakemake.input.raw[0])
clusters = pd.read_parquet(snakemake.input.cluster[0])
raw_cluster = raw.T.join(clusters).groupby('cluster').sum().T
raw_cluster.columns = [f'clus{x}' for x in raw_cluster.columns]
raw_cluster.to_parquet(snakemake.output[0])

import pandas as pd

tpm = pd.read_parquet(snakemake.input.tpm[0])
clusters = pd.read_parquet(snakemake.input.cluster[0])
tpm_cluster = tpm.T.join(clusters).groupby('cluster').sum().T
tpm_cluster.columns = [f'clus{x}' for x in tpm_cluster.columns]
tpm_cluster.to_parquet(snakemake.output[0])

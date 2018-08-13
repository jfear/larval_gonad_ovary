"""Aggregate individual counts matrices together."""
from pathlib import Path
import pandas as pd

dfs = []
for fname in snakemake.input:
    name = Path(fname).parents[0].name
    dat = pd.read_csv(fname, sep='\t', comment='#', usecols=[0, 6], index_col=0)
    dat.name = 'FBgn'
    dat.columns = [name]
    dfs.append(dat)

df = pd.concat(dfs, axis=1)
df.to_parquet(snakemake.output[0])

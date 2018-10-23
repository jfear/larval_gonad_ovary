import pandas as pd
from larval_gonad_ovary.normalization import zscore

df = pd.read_parquet(snakemake.input[0])
norm = zscore(df)
norm.to_parquet(snakemake.output[0])

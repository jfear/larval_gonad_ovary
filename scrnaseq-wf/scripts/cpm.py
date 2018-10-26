import pandas as pd
from larval_gonad_ovary.normalization import cpm

raw = pd.read_parquet(snakemake.input[0])
norm = cpm(raw)
norm.to_parquet(snakemake.output[0])

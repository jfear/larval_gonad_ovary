import pandas as pd
from larval_gonad_ovary.normalization import tpm

raw = pd.read_parquet(snakemake.input[0])
gene_lengths = pd.read_parquet(snakemake.input[1]).reindex(raw.index).gene_length
norm = tpm(raw, gene_lengths)
norm.to_parquet(snakemake.output[0])

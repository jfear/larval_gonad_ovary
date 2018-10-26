import pandas as pd


def pop(df, res):
    res = df[res]
    res.name = 'cluster'
    return res.to_frame()


df = pd.read_parquet(str(snakemake.input[0]))

res4 = pop(df, 'res.0.4')
res4.to_parquet(snakemake.output.res4[0])

res6 = pop(df, 'res.0.6')
res6.to_parquet(snakemake.output.res6[0])

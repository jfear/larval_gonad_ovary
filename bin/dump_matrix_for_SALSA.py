#!/usr/bin/env python
"""Quick script to dump read count matrices for use in SALSA.

This script needs to be run on a FAT node with ~100GB of RAM.
"""
from larval_gonad_ovary.logging import logger
from larval_gonad_ovary.io import cellranger_counts


def dump_counts(name):
    fname = f'../scrnaseq-wf/data/scrnaseq_samples/{name}/outs/raw_gene_bc_matrices_h5.h5'
    df = cellranger_counts(fname)

    #logger.info(f'Wirting {name} full matrix')
    #df.to_csv(f'../output/{name}_full_count_matrix_for_salsa.csv')

    logger.info(f'Wirting {name} barcode totals')
    df.sum(axis=0).to_csv(f'../output/{name}_barcode_count_matrix_for_salsa.csv')

    logger.info(f'Wirting {name} gene totals')
    df.sum(axis=1).to_csv(f'../output/{name}_gene_count_matrix_for_salsa.csv')


def main():
    for name in ['ovary1', 'ovary2', 'ovary3']:
        logger.info(f'Processing {name}')
        dump_counts(name)


if __name__ == '__main__':
    logger.info('Starting Script')
    main()
    logger.info('Script Complete')

"""Identification and selection of cells from 10x output.

The first processing step after running cell ranger is determining what
constitutes a cell and what is background. Cell ranger's algorithm does not
seem to work well in our situation due to differences in RNA content. The
number of UMI seems to be the best criteria to perform initial filtering,
perhaps followed by a gene filtering.

This module provides a set of functions for handling parsing of 10x HDF5
formats and munging data do decide what selection criteria should be.

"""
from collections import namedtuple

import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
import tables

from .config import memory

SOMA = [
]

EARLY_GERM = [
]

LATE_GERM = [
]

NUCS = ['A', 'C', 'G', 'T']
NUCS_INVERSE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

CellRangerCounts = namedtuple('CellRangerCounts',
                              ['matrix', 'gene_ids', 'barcodes'])


def compress_seq(s: str):
    """ Pack a DNA sequence (no Ns!) into a 2-bit format, in a 64-bit unit

    Most significant bit is set if there was an error

    Based on code from: https://github.com/10XGenomics/cellranger

    cellranger/lib/python/cellranger/utils.py

    """
    bits = 64
    assert len(s) <= (bits/2 - 1)
    result = 0
    for nuc in s:
        if nuc not in NUCS_INVERSE:
            return 1 << (bits - 1)
        result = result << 2
        result = result | NUCS_INVERSE[nuc]
    return result


def decompress_seq(x: int, length=16):
    """ Unpack a DNA sequence from a 2-bit format

    Based on code from: https://github.com/10XGenomics/cellranger

    cellranger/lib/python/cellranger/utils.py

    Parameters
    ----------
    x : int
        Number sequence to be decoded.
    length : int
        Length of the barcode. This can be found in the molecular info hdf5
        file from 10x genome.
        molInfo.get_node_attr('/metrics', 'chemistry_barcode_read_length')

    """
    bits = 64
    x = np.uint64(x)
    assert length <= (bits/2 - 1)
    if x & (1 << (bits-1)):
        return 'N' * length
    result = bytearray(length)
    for i in range(length):
        result[(length-1)-i] = bytearray(NUCS[x & np.uint64(0b11)].encode())[0]
        x = x >> np.uint64(2)
    return result.decode()


def two_bit_mapper(iterable):
    """Return a dictionary mapping 2bit encoded Seqs.

    Parameters
    ----------
    iterable : list-like
        Unique list of 2bit encoded sequences.

    Returns
    -------
    dict : Mapper from encoded to decoded

    """
    return {k: decompress_seq(k) for k in iterable}


def decode_cell_names(iterable):
    """Use two_bit_mapper to decode cell names.

    iterable : np.array
        An array of twobit encoded cell names.

    """
    mapper = two_bit_mapper(np.unique(iterable))
    return [mapper[x] for x in iterable]


@memory.cache
def cellranger_umi(fname):
    with tables.open_file(fname, 'r') as f:
        group = f.get_node('/')
        cell_ids = getattr(group, 'barcode').read()
        umi = getattr(group, 'umi').read()
        read_cnts = getattr(group, 'reads').read()

    cell_names = decode_cell_names(cell_ids)

    return pd.DataFrame(dict(
        cell_id=cell_names,
        umi=umi,
        read_cnt=read_cnts
    ))


@memory.cache
def cellranger_counts(fname, genome='dm6.16'):
    """Import cell ranger counts.

    Cell ranger stores it counts tables in a hdf5 formatted file. This reads
    this file and outputs them as a DataFrame.

    Parameters
    ----------
    fname : str
        Name of hdf5 store.
    genome : str
        Group where data is stored.
    barcodes : list of int
        Encoded barcodes names to filter by

    Returns
    -------
    namedtuple: matrix, gene_ids, barcodes

    """
    with tables.open_file(fname, 'r') as f:
        try:
            group = f.get_node(f.root, genome)
        except tables.NoSuchNodeError:
            print("That genome does not exist in this file.")
            return None
        gene_ids = getattr(group, 'genes').read()
        barcodes = getattr(group, 'barcodes').read()
        data = getattr(group, 'data').read()
        indices = getattr(group, 'indices').read()
        indptr = getattr(group, 'indptr').read()
        shape = getattr(group, 'shape').read()

    matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
    gene_ids = np.array([x.decode() for x in gene_ids])
    barcodes = np.array([x.decode().replace('-1', '') for x in barcodes])

    return CellRangerCounts(matrix, gene_ids, barcodes)


def calc_num_genes_on(cr: CellRangerCounts):
    """Number of genes that have >0 reads.

    Series containng the number of genes expressed by cell ID given a
    CellRangerCounts object.

    """
    num_genes_on = np.asarray((cr.matrix > 0).sum(axis=0))[0]
    idx = pd.Index(cr.barcodes, name='cell_id')
    return pd.Series(data=num_genes_on, index=idx, name='gene_cnt')


def cells_with_min_gene_expression(cr: CellRangerCounts, cutoff=200):
    """Get cell ids with more than `cutoff` genes.

    We are tyring to determine what kind of cell filtering should be preformed.
    In downstream analysis we use a gene level filter to remove cells that have
    fewer than `cutoff` expressed genes.  I am trying to determine if the cell
    ranger filter is having any affect on this or not.

    Parameters
    ----------
    cr : CellRangerCounts
        The raw_gene_bc_matrices_h5.h5 from cell ranger.
    cutoff : int
        The minimum number of expressed genes.

    Example
    -------
    >>> cr = cellranger_counts(...)
    >>> cells_on = cells_with_min_gene_expression(cr, cutoff=300)

    """
    num_genes_on = calc_num_genes_on(cr)
    cells_on = num_genes_on > cutoff
    return num_genes_on[cells_on].index.unique().tolist()


def filter_gene_counts_by_barcode(barcodes: np.array, cr: CellRangerCounts):
    """Given the raw gene counts matrix pull out specific cells.

    Example
    -------
    >>> cr = cellranger_counts(...)
    >>> barcodes = [...]
    >>> filter_gene_counts_by_barcode(barcodes, cr)

    """
    mask = np.in1d(cr.barcodes, barcodes)
    bcs = cr.barcodes[mask]
    matrix = cr.matrix[:, mask].todense()

    return pd.DataFrame(
        data=matrix,
        index=pd.Index(cr.gene_ids, name='FBgn'),
        columns=bcs
    )


def get_number_of_expressed_genes(fname):
    """Get number of genes with >0 reads.

    fname : str
        Path to the gene counts table.
    """
    cnts = pd.read_csv(fname, sep='\t')
    return (cnts.sum(axis=1) > 0).sum()


@memory.cache
def build_umi_gene_count_table(cr_raw, cr_umi):
    cr = cellranger_counts(cr_raw)
    umi = cellranger_umi(cr_umi)

    num_genes_on = calc_num_genes_on(cr)
    umi_cnts = umi.query('read_cnt > 0').groupby('cell_id').size().to_frame()
    umi_cnts.columns = ['umi_cnt']

    return umi_cnts.join(num_genes_on).sort_values('umi_cnt', ascending=False)

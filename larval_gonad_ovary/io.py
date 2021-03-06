"""Collection of io related items."""
import numpy as np
import pandas as pd
import scipy.sparse as sp_sparse
import tables

NUCS = ['A', 'C', 'G', 'T']
NUCS_INVERSE = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def compress_seq(s: str):
    """ Pack a DNA sequence (no Ns!) into a 2-bit format, in a 64-bit uint

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
    """ Un-pack a DNA sequence from a 2-bit format

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


def cellranger_umi(fname):
    with tables.open_file(fname, 'r') as f:
        group = f.get_node('/')
        cell_ids = getattr(group, 'barcode').read()
        umi = getattr(group, 'umi').read()
        read_cnts = getattr(group, 'reads').read()

    return pd.DataFrame(dict(
        cell_id=cell_ids,
        umi=umi,
        read_cnt=read_cnts
    ))


# Import Cell Ranger molecular info
def cellranger_counts(fname, genome='dm6.16', barcodes=None):
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

    """
    with tables.open_file(fname, 'r') as f:
        try:
            group = f.get_node(f.root, genome)
        except tables.NoSuchNodeError:
            print("That genome does not exist in this file.")
            return None
        gene_ids = getattr(group, 'genes').read()
        bcs = getattr(group, 'barcodes').read()
        data = getattr(group, 'data').read()
        indices = getattr(group, 'indices').read()
        indptr = getattr(group, 'indptr').read()
        shape = getattr(group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)

    if barcodes is not None:
        cell_ids = [bytes(decompress_seq(x) + '-1', 'utf-8') for x in barcodes]
        idx = []
        for i, bc in enumerate(bcs):
            if bc in cell_ids:
                idx.append(i)

        bcs = bcs[idx]
        matrix = matrix[:, idx]

    return pd.DataFrame(
        data=matrix.todense(),
        index=[x.decode() for x in gene_ids],
        columns=[x.decode() for x in bcs]
    )

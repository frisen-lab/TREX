from itertools import combinations
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
import pathlib
import re
import seaborn as sns
from tinyalign import hamming_distance


def load_reads(data_dir: pathlib.Path) -> pd.DataFrame:
    """Loads saved reads into a DataFrame"""
    READ_DIR = data_dir / 'reads.txt'
    return pd.read_csv(READ_DIR, delimiter='\t')


def load_molecules(data_dir: pathlib.Path) -> pd.DataFrame:
    """Loads saved molecules before correcting into a DataFrame."""
    MOLS_DIR = data_dir / 'molecules.txt'
    return pd.read_csv(MOLS_DIR, delimiter='\t')


def load_molecules_corrected(data_dir: pathlib.Path) -> pd.DataFrame:
    """Loads saved molecules after correcting into a DataFrame."""
    MOLS_DIR = data_dir / 'molecules_corrected.txt'
    return pd.read_csv(MOLS_DIR, delimiter='\t')


def load_cells(data_dir: pathlib.Path,
               filtered: bool = True) -> pd.DataFrame:
    cells_df = []

    filename = 'cells_filtered.txt' if filtered else 'cells.txt'

    with open(data_dir / filename) as file:
        next(file)
        while line := file.readline():
            cell_id, barcode_info = line[:-1].split('\t:\t')
            barcode_info = barcode_info.split('\t')

            barcodes_df = {'barcode': barcode_info[::2],
                           'counts': barcode_info[1::2]}
            barcodes_df = pd.DataFrame(barcodes_df,
                                       index=np.arange(len(barcode_info[::2])))
            barcodes_df['cell_id'] = cell_id
            cells_df.append(barcodes_df)

    cells_df = pd.concat(cells_df, ignore_index=True)
    cells_df['barcode'] = cells_df.barcode.astype('category')
    cells_df['counts'] = cells_df.counts.astype(int)
    cells_df['cell_id'] = cells_df.cell_id.astype('category')

    return cells_df


def load_umi_count_matrix(data_dir: pathlib.Path):
    """Loads saved UMI count matrix into a DataFrame."""
    UMI_DIR = data_dir / 'umi_count_matrix.csv'
    return pd.read_csv(UMI_DIR)


def load_clone_ids(data_dir: pathlib.Path):
    """Loads saved clone id and cell id into a DataFrame."""
    return pd.read_csv(data_dir / 'clones.txt')


def read_quality(reads: pd.DataFrame,
                 ax: plt.Axes = None,
                 add_description: bool = True) -> plt.Axes:
    """Plot histogram of how many times a molecule was read."""
    count_reads = reads.groupby(['#cell_id', 'umi']).agg('count')

    ax = sns.histplot(data=count_reads, x='clone_id', discrete=True, log=True,
                      ax=ax)
    plt.xlabel('Number of reads')
    plt.title('Number of reads per molecule')

    if add_description:
        txt = 'This plot shows how many times a molecule has been read. \n' \
              'The more reads per molecule, the better.'
        plt.text(0, -0.3, txt, transform=ax.transAxes, size=12)
    return ax


def get_read_per_molecule(df: pd.DataFrame) -> pd.Series:
    """Get a pandas Series with the number of reads per molecule."""
    return df.groupby(['#cell_id', 'umi']).clone_id.agg('count')


def get_length_read(df: pd.DataFrame,
                    molecules_dataframe: bool = True) -> pd.Series:
    """Get a pandas Series with the number of nucleotides read per molecule.
    molecules_dataframe is set to True by default, if a cells DataFrame is being
    used, then this must be False."""
    if molecules_dataframe:
        return df.clone_id.apply(lambda x: len(re.sub("[-0]", "", x)))
    else:
        return df.apply(lambda x: [len(
            re.sub("[-0]", "", x.barcode)), ] * x.counts,
                        axis=1).explode()


def get_barcodes_per_cell(df: pd.DataFrame,
                          molecules_dataframe: bool = True) -> pd.Series:
    """Get a pandas Series with the number of barcode molecules found per cell.
    molecules_dataframe is set to True by default, if a cells DataFrame is being
    used, then this must be False."""
    if molecules_dataframe:
        return df.groupby(['#cell_id']).umi.agg('count')
    else:
        return df.groupby(['cell_id']).counts.sum()


def get_molecules_per_barcodes(df: pd.DataFrame,
                               molecules_dataframe: bool = True) -> pd.Series:
    """Get a pandas Series with the number of barcode molecules found per unique
    barcode. molecules_dataframe is set to True by default, if a cells
    DataFrame is being used, then this must be False."""
    if molecules_dataframe:
        return df.groupby(['clone_id']).umi.agg('count')
    else:
        return df.groupby(['barcode']).counts.sum()


def get_unique_barcodes_per_cell(df: pd.DataFrame,
                                 molecules_dataframe: bool = True) -> pd.Series:
    """Get a pandas Series with the number of unique barcode molecules found per
    cell. molecules_dataframe is set to True by default, if a cells DataFrame is
    being used, then this must be False."""
    if molecules_dataframe:
        return df.groupby('#cell_id').clone_id.unique().apply(len)
    else:
        return df.groupby(['cell_id']).barcode.unique().apply(len)


def add_barcodes_per_clone(clones: pd.DataFrame,
                           cells_filtered: pd.DataFrame) -> pd.DataFrame:
    """Add  the unique barcode molecules found in each clone to the clones
    DataFrame. clones is the clones dataframe from clones.txt and cells_filtered
    is the dataframe containing barcodes per cell (cells_filtered.txt)."""
    cells_barcode = cells_filtered.groupby('cell_id').barcode.apply(set)
    clones['barcodes'] = clones.cell_id.apply(lambda x: cells_barcode[x])
    return clones


def get_barcodes_per_clone(clones: pd.DataFrame) -> pd.Series:
    """Get a pandas Series with the number of unique barcodes found per
    clone."""
    return clones.groupby('#clone_id').barcodes.apply(
        lambda x: len(set.union(*x)))


def get_clone_sizes(clones: pd.DataFrame) -> pd.Series:
    """Get a pandas Series with the umber of cells found in each
    clone. clones is the clones dataframe from clones.txt."""
    return clones.groupby('#clone_id').cell_id.count()


def plot_discrete_histogram(series: pd.Series,
                            title: str = None,
                            xlabel: str = None, ylabel: str = None,
                            txt:str = None,
                            ax: plt.Axes = None):
    """plots a discrete histogram with optional title, x and y labels, and added
     text below."""
    if ax is None:
        fig, ax = plt.subplots(1, 1)

    ax = sns.histplot(series, discrete=True, log=True, ax=ax)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    if txt is not None:
        plt.sca(ax)
        plt.text(0, -0.2, txt, transform=ax.transAxes, size=12, wrap=True,
                 ha='left', va='top')

    return ax


def length_read(molecules: pd.DataFrame,
                ax: plt.Axes = None,
                add_description: bool = True) -> plt.Axes:
    """Plot histogram of how many bases were adequately read per barcode."""
    molecules['stripped_barcode'] = molecules.clone_id.apply(
        lambda x: re.sub("[-0]", "", x))

    ax = sns.histplot(molecules.stripped_barcode.apply(len), discrete=True,
                      log=True, ax=ax)
    plt.xlabel('Number of detected bases')
    plt.title('Length of computed molecules')

    if add_description:
        txt = 'This plot shows how many bases have been read per molecule. \n' \
              'Ideally all 30 bases have been read. If very few bases are ' \
              'read, \nwe can not be sure how to complete the missing bases.'
        plt.text(0, -0.3, txt, transform=ax.transAxes, size=12)
    return ax


def molecules_per_cell(molecules: pd.DataFrame,
                       ax: plt.Axes = None,
                       add_description: bool = True) -> plt.Axes:
    """Plot histogram of how many molecules were detected per cell."""
    count_reads = molecules.groupby(['#cell_id']).umi.agg('count')

    ax = sns.histplot(count_reads.values, discrete=True, log=True, ax=ax)
    plt.xlabel('Molecules per cell')
    plt.title('Number of molecules')

    if add_description:
        txt = 'This plot shows how many molecules were found per cell. \n' \
              'Cell that have a single molecule might be filtered out later. \n' \
              'Cells that have too many might be result of non-removed \n doublets.'
        plt.text(0, -0.3, txt, transform=ax.transAxes, size=12)
    return ax


def molecules_per_barcode(molecules: pd.DataFrame,
                          ax: plt.Axes = None,
                          add_description : bool = True) -> plt.Axes:
    """Plot histogram of how many molecules were detected per viral barcode."""
    count_reads = molecules.groupby(['clone_id']).umi.agg('count')

    ax = sns.histplot(count_reads.values, discrete=True, log=True, ax=ax)
    plt.xlabel('Molecules per barcode')
    plt.title('Number of molecules')

    if add_description:
        txt = 'This plot shows how many molecules were found per barcode. \n' \
              'Barcodes that appear a few times might not be found in more \n' \
              'cells. Barcodes that have too many might be result of \n' \
              'contamination or alignment problems or big clones.'
        plt.text(0, -0.3, txt, transform=ax.transAxes, size=12)
    return ax


def unique_barcodes_per_cell(molecules: pd.DataFrame,
                             ax: plt.Axes = None,
                             add_description : bool = True) -> plt.Axes:
    """Plot histogram of how many unique barcodes were detected per cell."""
    count_reads = molecules.groupby('#cell_id').clone_id.unique().apply(len)

    ax = sns.histplot(count_reads.values, discrete=True, log=True, ax=ax)
    plt.xlabel('Unique barcodes per cell')
    plt.title('Number of unique barcodes')

    if add_description:
        txt = 'This plot shows how many unique barcodes were detected per cell.\n' \
              'Cells with many unique barcodes show either lots of infection \n' \
              'events or possible unfiltered doublets.'
        plt.text(0, -0.3, txt, transform=ax.transAxes, size=12)
    return ax


def hamming_distance_histogram(molecules: pd.DataFrame,
                               ax: plt.Axes = None,
                               ignore_incomplete: bool = True) -> plt.Axes:
    """Plot histogram of Hamming distance between barcodes. ignore_incomplete is
     set to True by default and it removes incomplete barcodes."""
    if ignore_incomplete:
        molecules = molecules[~molecules.clone_id.str.contains('-|0')]
    this_clone_ids = molecules.clone_id.unique()
    hamming_distances = np.empty([len(this_clone_ids)] * 2)

    def my_iter(barcode_list):
        for inds in combinations(np.arange(len(barcode_list)), 2):
            yield inds, barcode_list[inds[0]], barcode_list[inds[1]]

    # Hamming distance function
    def is_similar(args):
        inds, s, t = args
        if not ignore_incomplete:
            bad_chars = {'-', '0'}
            if bad_chars & set(s) or bad_chars & set(t):
                # Remove suffix and/or prefix where sequences do not overlap
                s = s.lstrip("-0")
                t = t[-len(s):]
                t = t.lstrip("-0")
                s = s[-len(t):]
                s = s.rstrip("-0")
                t = t[: len(s)]
                t = t.rstrip("-0")
                s = s[: len(t)]

        return inds, hamming_distance(s, t)

    for ind, val in map(is_similar, my_iter(this_clone_ids)):
        hamming_distances[ind[0], ind[1]] = val

    vals = hamming_distances[np.triu_indices_from(hamming_distances, 1)]
    ax = sns.histplot(vals, discrete=True, log=True, ax=ax)

    ax.set_title('Hamming Distance Histogram')
    ax.set_xlabel('Hamming Distance')

    return ax


def jaccard(cell_1: npt.ArrayLike, cell_2: npt.ArrayLike) -> float:
    """Estimates the Jaccard similarity between two boolean vectors taken from
    the UMI count matrix."""
    return sum(np.logical_and(cell_1, cell_2)) / sum(np.logical_or(cell_1,
                                                                   cell_2))


def jaccard_similarity_matrix(umi_count: pd.DataFrame) -> npt.ArrayLike:
    """Builds Jaccard similarity matrix from the UMI count matrix loaded as
    pandas DataFrame.

    Note: columns are cell IDs and first column is disregarded as it usually has
    the index to barcodes"""
    this_cell_ids = umi_count.columns[1:]
    jaccard_matrix = np.zeros([len(this_cell_ids)] * 2)

    # Iterator over combinations of cells
    def my_iter(cell_id_list):
        for inds in combinations(np.arange(len(cell_id_list)), 2):
            clones_cell_1 = umi_count[cell_id_list[inds[0]]].values > 0
            clones_cell_2 = umi_count[cell_id_list[inds[1]]].values > 0
            yield inds, clones_cell_1, clones_cell_2

    def my_jaccard(args):
        inds, clones_cell_1, clones_cell_2 = args
        return inds, jaccard(clones_cell_1, clones_cell_2)

    for ind, val in map(my_jaccard, my_iter(this_cell_ids)):
        jaccard_matrix[ind[0], ind[1]] = val

    jaccard_matrix = jaccard_matrix + jaccard_matrix.T
    jaccard_matrix[np.diag_indices_from(jaccard_matrix)] = 1

    return jaccard_matrix


def jaccard_histogram(jaccard_matrix: npt.ArrayLike,
                      ax: plt.Axes = None) -> plt.Axes:
    """Plots the Jaccard similarity histogram between cells."""
    vals = jaccard_matrix[np.triu_indices_from(jaccard_matrix, 1)]
    ax = sns.histplot(vals, log=True, ax=ax)

    ax.set_title('Jaccard Similarity between cells Histogram')
    ax.set_xlabel('Jaccard Similarity')

    return ax


def plot_jaccard_matrix(jaccard_matrix: npt.ArrayLike,
                        ax: plt.Axes = None) -> plt.Axes:
    """Orders with reverse Cuthill-McKee algorithm the cells in the Jaccard
     Similarity matrix and plots it for visualization."""
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import reverse_cuthill_mckee

    graph = csr_matrix(jaccard_matrix)
    swaps = reverse_cuthill_mckee(graph, symmetric_mode=True)

    jaccard_matrix = jaccard_matrix[swaps]
    jaccard_matrix = jaccard_matrix[:, swaps]

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(11, 10))
    im = ax.imshow(jaccard_matrix, interpolation='none', cmap='Blues',
                   rasterized=True, vmin=0, vmax=1)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylabel('Cell ID')
    ax.set_xlabel('Cell ID')
    plt.colorbar(im, ax=ax)

    return ax

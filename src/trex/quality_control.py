import pathlib
import re
from itertools import combinations

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
import seaborn as sns
from tinyalign import hamming_distance


def read_cells(data_dir: pathlib.Path, filtered: bool = True) -> pd.DataFrame:
    cells_df = []

    filename = "cells_filtered.txt" if filtered else "cells.txt"

    with open(data_dir / filename) as file:
        next(file)
        for line in file:
            cell_id, clone_id_info = line[:-1].split("\t:\t")
            clone_id_info = clone_id_info.split("\t")

            clone_ids_df = {
                "clone_id": clone_id_info[::2],
                "counts": clone_id_info[1::2],
            }
            clone_ids_df = pd.DataFrame(
                clone_ids_df, index=np.arange(len(clone_id_info[::2]))
            )
            clone_ids_df["cell_id"] = cell_id
            cells_df.append(clone_ids_df)

    cells_df = pd.concat(cells_df, ignore_index=True)
    cells_df["clone_id"] = cells_df.clone_id.astype("category")
    cells_df["counts"] = cells_df.counts.astype(int)
    cells_df["cell_id"] = cells_df.cell_id.astype("category")

    return cells_df


def read_umi_count_matrix(data_dir: pathlib.Path) -> pd.DataFrame:
    """Loads saved UMI count matrix into a DataFrame."""
    umi_dir = data_dir / "umi_count_matrix.csv"
    umi_matrix = pd.read_csv(umi_dir)
    umi_matrix.set_index("Unnamed: 0", inplace=True, verify_integrity=True)
    return umi_matrix


def read_clone_ids(data_dir: pathlib.Path):
    """Loads saved clone id and cell id into a DataFrame."""
    df = pd.read_table(data_dir / "clones.txt")

    # In case previous version of clones is loaded, clone_nr was called #clone_id
    if "clone_nr" not in df.columns:
        df.rename(columns={"#clone_id": "clone_nr"}, inplace=True)
    return df


def plot_read_quality(
    reads: pd.DataFrame, ax: plt.Axes = None, add_description: bool = True
) -> plt.Axes:
    """Plot histogram of how many times a molecule was read."""
    count_reads = reads.groupby(["#cell_id", "umi"]).agg("count")

    ax = sns.histplot(data=count_reads, x="clone_id", discrete=True, log=True, ax=ax)
    plt.xlabel("Number of reads")
    plt.title("Number of reads per molecule")

    if add_description:
        txt = (
            "This plot shows how many times a molecule has been read. \n"
            "The more reads per molecule, the better."
        )
        plt.text(0, -0.3, txt, transform=ax.transAxes, size=12)
    return ax


def get_read_per_molecule(df: pd.DataFrame) -> pd.Series:
    """Get a pandas Series with the number of reads per molecule."""
    return df.groupby(["#cell_id", "umi"]).clone_id.agg("count")


def get_nucleotides_per_molecule(df: pd.DataFrame) -> pd.Series:
    """Get a pandas Series with the number of nucleotides read per molecule."""
    if "counts" in df.columns:
        return df.apply(
            lambda x: [
                len(re.sub("[-0]", "", x.clone_id)),
            ]
            * x.counts,
            axis=1,
        ).explode()
    else:
        return df.clone_id.apply(lambda x: len(re.sub("[-0]", "", x)))


def get_clone_ids_per_cell(
    df: pd.DataFrame, molecules_dataframe: bool = True
) -> pd.Series:
    """
    Get a pandas Series with the number of cloneID molecules found per cell.
    molecules_dataframe is set to True by default, if a cells DataFrame is being
    used, then this must be False.
    """
    if molecules_dataframe:
        return df.groupby(["#cell_id"]).umi.agg("count")
    else:
        return df.groupby(["cell_id"]).counts.sum()


def get_molecules_per_clone_ids(
    df: pd.DataFrame, molecules_dataframe: bool = True
) -> pd.Series:
    """
    Get a pandas Series with the number of cloneID molecules found per unique
    cloneID. molecules_dataframe is set to True by default, if a cells
    DataFrame is being used, then this must be False.
    """
    if molecules_dataframe:
        return df.groupby(["clone_id"]).umi.agg("count")
    else:
        return df.groupby(["clone_id"]).counts.sum()


def get_unique_clone_ids_per_cell(
    df: pd.DataFrame, molecules_dataframe: bool = True
) -> pd.Series:
    """
    Get a pandas Series with the number of unique cloneID molecules found per
    cell. molecules_dataframe is set to True by default, if a cells DataFrame is
    being used, then this must be False.
    """
    if molecules_dataframe:
        return df.groupby("#cell_id").clone_id.unique().apply(len)
    else:
        return df.groupby(["cell_id"]).clone_id.unique().apply(len)


def add_clone_ids_per_clone(
    clones: pd.DataFrame, cells_filtered: pd.DataFrame
) -> pd.DataFrame:
    """
    Add  the unique cloneID molecules found in each clone to the clones
    DataFrame. clones is the clones dataframe from clones.txt and cells_filtered
    is the dataframe containing cloneIDs per cell (cells_filtered.txt).
    """
    cells_clone_id = cells_filtered.groupby("cell_id").clone_id.apply(set)
    clones["clone_ids"] = clones.cell_id.apply(lambda x: cells_clone_id[x])
    return clones


def get_clone_ids_per_clone(clones: pd.DataFrame) -> pd.Series:
    """
    Get a pandas Series with the number of unique cloneIDs found per
    clone.
    """
    return clones.groupby("clone_nr").clone_ids.apply(lambda x: len(set.union(*x)))


def get_clone_sizes(clones: pd.DataFrame) -> pd.Series:
    """
    Get a pandas Series with the number of cells found in each
    clone. clones is the clones dataframe from clones.txt.
    """
    return clones.groupby("clone_nr").cell_id.count()


def plot_discrete_histogram(
    series: pd.Series,
    title: str = None,
    xlabel: str = None,
    ylabel: str = None,
    text: str = None,
    axes: plt.Axes = None,
):
    """
    Plots a discrete histogram with optional title, x and y labels, and added
    text below.
    """
    if axes is None:
        fig, axes = plt.subplots(1, 1)

    axes = sns.histplot(series, discrete=True, log=True, ax=axes)

    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    axes.set_title(title)

    if text is not None:
        plt.sca(axes)
        plt.text(
            0,
            -0.2,
            text,
            transform=axes.transAxes,
            size=12,
            wrap=True,
            ha="left",
            va="top",
        )

    return axes


def plot_molecules_per_cell(
    molecules: pd.DataFrame, ax: plt.Axes = None, add_description: bool = True
) -> plt.Axes:
    """Plot histogram of how many molecules were detected per cell."""
    count_reads = molecules.groupby(["#cell_id"]).umi.agg("count")

    ax = sns.histplot(count_reads.values, discrete=True, log=True, ax=ax)
    plt.xlabel("Molecules per cell")
    plt.title("Number of molecules")

    if add_description:
        txt = (
            "This plot shows how many molecules were found per cell. \n"
            "Cell that have a single molecule might be filtered out later. \n"
            "Cells that have too many might be result of non-removed \n doublets."
        )
        plt.text(0, -0.3, txt, transform=ax.transAxes, size=12)
    return ax


def plot_molecules_per_clone_id(
    molecules: pd.DataFrame, ax: plt.Axes = None, add_description: bool = True
) -> plt.Axes:
    """Plot histogram of how many molecules were detected per viral cloneID."""
    count_reads = molecules.groupby(["clone_id"]).umi.agg("count")

    ax = sns.histplot(count_reads.values, discrete=True, log=True, ax=ax)
    plt.xlabel("Molecules per cloneID")
    plt.title("Number of molecules")

    if add_description:
        txt = (
            "This plot shows how many molecules were found per cloneID. \n"
            "cloneIDs that appear a few times might not be found in more \n"
            "cells. cloneIDs that have too many might be result of \n"
            "contamination or alignment problems or big clones."
        )
        plt.text(0, -0.3, txt, transform=ax.transAxes, size=12)
    return ax


def plot_unique_clone_ids_per_cell(
    molecules: pd.DataFrame, ax: plt.Axes = None, add_description: bool = True
) -> plt.Axes:
    """Plot histogram of how many unique cloneIDs were detected per cell."""
    count_reads = molecules.groupby("#cell_id").clone_id.unique().apply(len)

    ax = sns.histplot(count_reads.values, discrete=True, log=True, ax=ax)
    plt.xlabel("Unique cloneIDs per cell")
    plt.title("Number of unique cloneIDs")

    if add_description:
        txt = (
            "This plot shows how many unique cloneIDs were detected per cell.\n"
            "Cells with many unique cloneIDs show either lots of infection \n"
            "events or possible unfiltered doublets."
        )
        plt.text(0, -0.3, txt, transform=ax.transAxes, size=12)
    return ax


def plot_hamming_distance_histogram(
    molecules: pd.DataFrame, ax: plt.Axes = None, ignore_incomplete: bool = True
) -> plt.Axes:
    """
    Plot histogram of Hamming distance between cloneIDs. ignore_incomplete is
    set to True by default and it removes incomplete cloneIDs.
    """
    if ignore_incomplete:
        molecules = molecules[~molecules.clone_id.str.contains("-|0")]
    this_clone_ids = molecules.clone_id.unique()
    hamming_distances = np.empty([len(this_clone_ids)] * 2)

    def my_iter(clone_id_list):
        for inds in combinations(np.arange(len(clone_id_list)), 2):
            yield inds, clone_id_list[inds[0]], clone_id_list[inds[1]]

    # Hamming distance function
    def is_similar(args):
        inds, s, t = args
        if not ignore_incomplete:
            bad_chars = {"-", "0"}
            if bad_chars & set(s) or bad_chars & set(t):
                # Remove suffix and/or prefix where sequences do not overlap
                s = s.lstrip("-0")
                t = t[-len(s) :]
                t = t.lstrip("-0")
                s = s[-len(t) :]
                s = s.rstrip("-0")
                t = t[: len(s)]
                t = t.rstrip("-0")
                s = s[: len(t)]

        return inds, hamming_distance(s, t)

    for ind, val in map(is_similar, my_iter(this_clone_ids)):
        hamming_distances[ind[0], ind[1]] = val

    vals = hamming_distances[np.triu_indices_from(hamming_distances, 1)]
    ax = sns.histplot(vals, discrete=True, log=True, ax=ax)

    ax.set_title("Hamming Distance Histogram")
    ax.set_xlabel("Hamming Distance")

    return ax


def jaccard(cell_1: npt.ArrayLike, cell_2: npt.ArrayLike) -> float:
    """
    Estimate the Jaccard similarity between two boolean vectors taken from
    the UMI count matrix.
    """
    return sum(np.logical_and(cell_1, cell_2)) / sum(np.logical_or(cell_1, cell_2))


def jaccard_similarity_matrix(umi_count: pd.DataFrame) -> npt.ArrayLike:
    """
    Build Jaccard similarity matrix from the UMI count matrix loaded as
    pandas DataFrame.

    Note: columns are cell IDs and first column is disregarded as it usually has
    the index to cloneIDs
    """
    this_cell_ids = umi_count.columns[1:]
    jaccard_matrix = np.zeros([len(this_cell_ids)] * 2)
    bool_umi = pd.DataFrame(
        umi_count.values > 0, columns=umi_count.columns, index=umi_count.index
    )
    cell_id_position_dict = {
        this_cell_id: i for i, this_cell_id in enumerate(this_cell_ids)
    }

    for id1, id2 in combinations(this_cell_ids, 2):
        clones_cell_1 = bool_umi[id1].values
        clones_cell_2 = bool_umi[id2].values
        jaccard_matrix[
            cell_id_position_dict[id1], cell_id_position_dict[id2]
        ] = jaccard(clones_cell_1, clones_cell_2)

    jaccard_matrix = jaccard_matrix + jaccard_matrix.T
    jaccard_matrix[np.diag_indices_from(jaccard_matrix)] = 1

    return jaccard_matrix


def jaccard_histogram(jaccard_matrix: npt.ArrayLike, axes: plt.Axes = None) -> plt.Axes:
    """Plot the Jaccard similarity histogram between cells."""
    vals = jaccard_matrix[np.triu_indices_from(jaccard_matrix, 1)]
    axes = sns.histplot(vals, log=True, ax=axes)

    axes.set_title("Jaccard Similarity between cells Histogram")
    axes.set_xlabel("Jaccard Similarity")

    return axes


def plot_jaccard_matrix(
    jaccard_matrix: npt.ArrayLike, axes: plt.Axes = None
) -> plt.Axes:
    """
    Order with reverse Cuthill-McKee algorithm the cells in the Jaccard
    Similarity matrix and plots it for visualization.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import reverse_cuthill_mckee

    graph = csr_matrix(jaccard_matrix)
    swaps = reverse_cuthill_mckee(graph, symmetric_mode=True)

    jaccard_matrix = jaccard_matrix[swaps]
    jaccard_matrix = jaccard_matrix[:, swaps]

    if axes is None:
        fig, axes = plt.subplots(1, 1, figsize=(11, 10))
    im = axes.imshow(
        jaccard_matrix,
        interpolation="none",
        cmap="Blues",
        rasterized=True,
        vmin=0,
        vmax=1,
    )

    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_ylabel("Cell ID")
    axes.set_xlabel("Cell ID")
    plt.colorbar(im, ax=axes)

    return axes

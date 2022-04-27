from pathlib import Path
from typing import List
from .cell import Cell
import operator
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings("ignore", "Conversion of the second argument of issubdtype")
    import loompy
import numpy as np


def write_count_matrix(path: Path, cells: List[Cell]):
    """Create a Read-count matrix with cells as columns and cloneIDs as rows"""
    clone_ids = set()
    for cell in cells:
        clone_ids.update(clone_id for clone_id in cell.counts)
    clone_ids = sorted(clone_ids)
    all_counts = [cell.counts for cell in cells]
    with open(path, "w") as f:
        f.write(",")
        f.write(",".join(cell.cell_id for cell in cells))
        f.write("\n")
        for clone_id in clone_ids:
            f.write(clone_id)
            f.write(",")
            values = [lic.get(clone_id, 0) for lic in all_counts]
            f.write(",".join(str(v) for v in values))
            f.write("\n")


def write_cells(path: Path, cells: List[Cell]) -> None:
    """Write cells to a tab-separated file"""
    with open(path, "w") as f:
        print(
            "#cell_id",
            ":",
            "clone_id1",
            "count1",
            "clone_id2",
            "count2",
            "...",
            sep="\t",
            file=f,
        )
        for cell in cells:
            row = [cell.cell_id, ":"]
            sorted_clone_ids = sorted(
                cell.counts, key=lambda x: cell.counts[x], reverse=True
            )
            if not sorted_clone_ids:
                continue
            for clone_id in sorted_clone_ids:
                row.extend([clone_id, cell.counts[clone_id]])
            print(*row, sep="\t", file=f)


def write_reads_or_molecules(path, mols_or_reads, require_umis=True, sort=True):
    with open(path, "w") as f:
        if require_umis:
            if sort:
                mols_or_reads = sorted(
                    mols_or_reads,
                    key=lambda mol_or_read: (
                        mol_or_read.umi,
                        mol_or_read.cell_id,
                        mol_or_read.clone_id,
                    ),
                )
            print("#cell_id", "umi", "clone_id", sep="\t", file=f)
            for mol_or_read in mols_or_reads:
                print(
                    mol_or_read.cell_id,
                    mol_or_read.umi,
                    mol_or_read.clone_id,
                    sep="\t",
                    file=f,
                )
        else:
            if sort:
                mols_or_reads = sorted(
                    mols_or_reads,
                    key=lambda mol_or_read: (mol_or_read.clone_id, mol_or_read.cell_id),
                )
            print("#cell_id", "clone_id", sep="\t", file=f)
            for mol_or_read in mols_or_reads:
                print(mol_or_read.cell_id, mol_or_read.clone_id, sep="\t", file=f)


def write_loom(cells: List[Cell], cellranger, output_dir, clone_id_length, top_n=6):
    """
    Create a loom file from a Cell Ranger result directory and augment it with information about
    the most abundant cloneIDs and their counts.
    """
    # For each cell, collect the most abundant cloneIDs and their counts
    # Maps cell_id to a list of (clone_id, count) pairs that represent the most abundant cloneIDs.
    most_abundant = dict()
    for cell in cells:
        if not cell.counts:
            continue
        counts = sorted(cell.counts.items(), key=operator.itemgetter(1))
        counts.reverse()
        counts = counts[:top_n]
        most_abundant[cell.cell_id] = counts

    loompy.create_from_cellranger(cellranger.sample_dir, outdir=output_dir)
    # create_from_cellranger() does not tell us the name of the created file,
    # so we need to re-derive it from the sample name.
    sample_name = cellranger.sample_dir.name
    loom_path = output_dir / (sample_name + ".loom")

    with loompy.connect(loom_path) as ds:
        # Cell ids in the loom file are prefixed by the sample name and a ':'. Remove that prefix.
        loom_cell_ids = [cell_id[len(sample_name) + 1 :] for cell_id in ds.ca.CellID]

        # Transform cloneIDs and count data
        # brings cloneID data into correct format for loom file.
        # Array must have same shape as all_cellIDs
        clone_id_lists = [[] for _ in range(top_n)]
        count_lists = [[] for _ in range(top_n)]
        for cell_id in loom_cell_ids:
            clone_id_counts = most_abundant.get(cell_id, [])
            # Fill up to a constant length
            while len(clone_id_counts) < top_n:
                clone_id_counts.append(("-", 0))

            for i, (clone_id, count) in enumerate(clone_id_counts):
                clone_id_lists[i].append(clone_id)
                count_lists[i].append(count)

        # Add cloneID and count information to loom file
        for i in range(top_n):
            ds.ca[f"cloneid_{i+1}"] = np.array(
                clone_id_lists[i], dtype="S%r" % clone_id_length
            )
            ds.ca[f"cloneid_count_{i+1}"] = np.array(count_lists[i], dtype=int)

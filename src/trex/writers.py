from pathlib import Path
from typing import List
from .cell import Cell


def write_count_matrix(path: Path, cells: List[Cell]):
    """Create a Read-count matrix with cells as columns and clone IDs as rows"""
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
    with open(path, 'w') as f:
        print(
            "#cell_id", ":", "clone_id1", "count1", "clone_id2", "count2", "...", sep="\t", file=f)
        for cell in cells:
            row = [cell.cell_id, ':']
            sorted_clone_ids = sorted(
                cell.counts, key=lambda x: cell.counts[x], reverse=True)
            if not sorted_clone_ids:
                continue
            for clone_id in sorted_clone_ids:
                row.extend([clone_id, cell.counts[clone_id]])
            print(*row, sep='\t', file=f)


def write_reads(path, reads, require_umis=True):
    with open(path, 'w') as f:
        if require_umis:
            print("#cell_id", "umi", "clone_id", sep="\t", file=f)
            for read in sorted(reads, key=lambda read: (read.umi, read.cell_id, read.clone_id)):
                print(read.cell_id, read.umi, read.clone_id, sep='\t', file=f)
        else:
            print("#cell_id", "clone_id", sep="\t", file=f)
            for read in sorted(reads, key=lambda read: (read.clone_id, read.cell_id)):
                print(read.cell_id, read.clone_id, sep='\t', file=f)


def write_molecules(path, molecules):
    with open(path, 'w') as f:
        print("#cell_id", "umi", "clone_id", sep="\t", file=f)
        for molecule in molecules:
            print(molecule.cell_id, molecule.umi, molecule.clone_id, sep='\t', file=f)
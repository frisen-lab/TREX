"""
Run on single cell 10X Chromium or spatial Visium data processed by Cell / Space Ranger software
"""
import re
import sys
import logging
import dataclasses
from pathlib import Path
from collections import Counter, defaultdict
from typing import List, Dict, Iterable, Optional, DefaultDict, Union

import pandas as pd

from . import (
    CommandLineError,
    add_file_logging,
    make_output_dir,
    add_common_arguments,
)
from .. import __version__
from ..cellranger import make_cellranger, CellRangerError
from ..writers import (
    write_count_matrix,
    write_cells,
    write_reads_or_molecules,
    write_loom,
)
from ..clustering import cluster_sequences
from ..clone import CloneGraph
from ..filters import is_low_complexity
from ..molecule import Molecule, compute_molecules
from ..cell import Cell, compute_cells
from ..error import TrexError
from ..dataset import DatasetReader


__author__ = "leonie.von.berlin@ki.se"

logger = logging.getLogger(__name__)


def main(args):
    output_dir = args.output
    try:
        make_output_dir(output_dir, args.delete)
    except FileExistsError:
        raise CommandLineError(
            f"Output directory '{output_dir}' already exists "
            "(use --delete to force deleting an existing output directory)"
        )

    add_file_logging(output_dir / "log.txt")
    logger.info(f"Trex {__version__}")
    logger.info("Command line arguments: %s", " ".join(sys.argv[1:]))

    allowed_cell_ids = None
    if args.filter_cellids:
        allowed_cell_ids = read_allowed_cellids(args.filter_cellids)
    excluded_clone_ids = None
    if args.filter_cloneids:
        excluded_clone_ids = read_excluded_clone_ids(args.filter_cloneids)
    transcriptome_inputs = args.path
    if args.samples:
        sample_names = args.samples.split(",")
    elif len(transcriptome_inputs) == 1:
        sample_names = [None]  # Do not modify suffixes
    else:
        sample_names = [path.name for path in transcriptome_inputs]
        logger.info("Using these sample names: %s", ", ".join(sample_names))
    if len(sample_names) != len(transcriptome_inputs):
        raise CommandLineError(
            "The number of sample names (--samples) must match the number of "
            "provided transcriptome datasets"
        )
    if args.amplicon:
        amplicon_inputs = args.amplicon
        if len(transcriptome_inputs) != len(amplicon_inputs):
            raise CommandLineError(
                "As many amplicon as transcriptome datasets must be provided"
            )
    else:
        amplicon_inputs = []

    highlight_cell_ids = []
    if args.highlight:
        with open(args.highlight) as f:
            highlight_cell_ids = [line.strip() for line in f]

    try:
        run_trex(
            output_dir,
            genome_name=args.genome_name,
            allowed_cell_ids=allowed_cell_ids,
            excluded_clone_ids=excluded_clone_ids,
            chromosome=args.chromosome,
            start=args.start - 1 if args.start is not None else None,
            end=args.end,
            transcriptome_inputs=transcriptome_inputs,
            amplicon_inputs=amplicon_inputs,
            sample_names=sample_names,
            prefix=args.prefix,
            correct_per_cell=args.correct_per_cell,
            max_hamming=args.max_hamming,
            min_length=args.min_length,
            jaccard_threshold=args.jaccard_threshold,
            keep_single_reads=args.keep_single_reads,
            keep_doublets=args.keep_doublets,
            should_write_umi_matrix=args.umi_matrix,
            should_run_visium=args.visium,
            should_plot=args.plot,
            highlight_cell_ids=highlight_cell_ids,
            should_write_loom=args.loom,
        )
    except (CellRangerError, TrexError) as e:
        raise CommandLineError(e)


def add_arguments(parser):
    groups = add_common_arguments(parser, smartseq=False)

    groups.filter.add_argument(
        "--keep-single-reads",
        action="store_true",
        default=False,
        help="Keep cloneIDs supported by only a single read. Default: Discard them",
    )
    groups.filter.add_argument(
        "--keep-doublets",
        action="store_true",
        default=False,
        help="Keep doublets. Default: Detect and remove doublets",
    )
    groups.filter.add_argument(
        "--visium",
        default=False,
        action="store_true",
        help="Adjust filter settings for 10x Visium data: Filter out cloneIDs only based on "
        "one read, but keep those with only one UMI",
    )
    groups.filter.add_argument(
        "--per-cell",
        help="Use only cloneIDs within the same cell for cloneID correction. "
        "Default: Use cloneIDs from all cells",
        default=False,
        action="store_true",
        dest="correct_per_cell",
    )

    groups.output.add_argument(
        "-l",
        "--loom",
        help="Create also a loom-file from Cell Ranger and clone data. "
        "File will have the same name as the run. Default: do not create a loom file",
        action="store_true",
    )
    groups.output.add_argument(
        "--umi-matrix",
        default=False,
        action="store_true",
        help="Create a UMI count matrix 'umi_count_matrix.csv' with "
        "cells as columns and cloneIDs as rows",
    )


def run_trex(
    output_dir: Path,
    *,
    transcriptome_inputs: List[Union[Path, str]],
    amplicon_inputs: List[Path],
    genome_name: Optional[str] = None,
    allowed_cell_ids: Optional[List[str]] = None,
    excluded_clone_ids: Optional[List[str]] = None,
    chromosome: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    sample_names: Optional[List[str]] = None,
    prefix: bool = False,
    max_hamming: int = 5,
    correct_per_cell: bool = False,
    min_length: int = 20,
    jaccard_threshold: float = 0.7,
    keep_single_reads: bool = False,
    keep_doublets: bool = False,
    should_write_umi_matrix: bool = False,
    should_run_visium: bool = False,
    should_plot: bool = False,
    highlight_cell_ids: Optional[List[str]] = None,
    should_write_loom: bool = False,
):
    if sample_names is not None and len(sample_names) != len(set(sample_names)):
        raise TrexError("The sample names need to be unique")

    dataset_reader = DatasetReader(
        output_dir, genome_name, chromosome, start, end, prefix
    )
    reads = dataset_reader.read_all(
        transcriptome_inputs, amplicon_inputs, sample_names, allowed_cell_ids
    )
    if not reads:
        raise TrexError("No reads left after --filter-cellids filtering")

    clone_ids = [
        r.clone_id for r in reads if "-" not in r.clone_id and "0" not in r.clone_id
    ]
    logger.info(
        f"Read {len(reads)} reads containing (parts of) the cloneID "
        f"({len(clone_ids)} full cloneIDs, {len(set(clone_ids))} unique)"
    )

    write_reads_or_molecules(output_dir / "reads.txt", reads)

    molecules = compute_molecules(reads)
    clone_ids = [
        m.clone_id for m in molecules if "-" not in m.clone_id and "0" not in m.clone_id
    ]
    logger.info(
        f"Detected {len(molecules)} molecules ({len(clone_ids)} full cloneIDs, "
        f"{len(set(clone_ids))} unique)"
    )
    write_reads_or_molecules(output_dir / "molecules.txt", molecules, sort=False)

    molecules = [m for m in molecules if not is_low_complexity(m.clone_id)]
    logger.info(f"{len(molecules)} remain after low-complexity filtering")
    write_reads_or_molecules(
        output_dir / "molecules_filtered.txt", molecules, sort=False
    )

    if correct_per_cell:
        corrected_molecules = correct_clone_ids_per_cell(
            molecules, max_hamming, min_length
        )
    else:
        corrected_molecules = correct_clone_ids(molecules, max_hamming, min_length)
    clone_ids = [
        m.clone_id
        for m in corrected_molecules
        if "-" not in m.clone_id and "0" not in m.clone_id
    ]
    logger.info(
        f"After cloneID correction, {len(set(clone_ids))} unique cloneIDs remain"
    )

    write_reads_or_molecules(
        output_dir / "molecules_corrected.txt", corrected_molecules, sort=False
    )

    if excluded_clone_ids is not None:
        similarity_set = SimilaritySet(excluded_clone_ids)
        corrected_molecules = [
            molecule
            for molecule in corrected_molecules
            if not similarity_set.contains(molecule.clone_id)
        ]
        clone_ids = [
            m.clone_id
            for m in corrected_molecules
            if "-" not in m.clone_id and "0" not in m.clone_id
        ]
        logger.info(
            f"After filtering cloneIDs, {len(set(clone_ids))} unique cloneIDs remain"
        )

    cells = compute_cells(corrected_molecules, min_length)
    logger.info(f"Detected {len(cells)} cells")
    write_cells(output_dir / "cells.txt", cells)

    if should_run_visium:
        cells = filter_visium(cells, corrected_molecules)
        logger.info(f"{len(cells)} filtered cells remain")
        write_cells(output_dir / "cells_filtered.txt", cells)
    else:
        cells = filter_cells(cells, corrected_molecules, keep_single_reads)
        logger.info(f"{len(cells)} filtered cells remain")
        write_cells(output_dir / "cells_filtered.txt", cells)

    if should_write_umi_matrix:
        logger.info("Writing UMI matrix")
        write_count_matrix(output_dir / "umi_count_matrix.csv", cells)

    clone_graph = CloneGraph(cells, jaccard_threshold=jaccard_threshold)

    with open(output_dir / "components.txt", "w") as components_file:
        print(
            clone_graph.components_txt(highlight_cell_ids), file=components_file, end=""
        )

    bridges = clone_graph.bridges()
    logger.info(f"Removing {len(bridges)} bridges from the graph")
    clone_graph.remove_edges(bridges)

    if not keep_doublets:
        doublets = clone_graph.doublets()
        logger.info(f"Removing {len(doublets)} doublets from the graph (first round)")
        clone_graph.remove_nodes(doublets)

        bridges2 = clone_graph.bridges()
        logger.info(f"Removing {len(bridges2)} bridges from the graph (second round)")
        clone_graph.remove_edges(bridges2)
        doublets2 = clone_graph.doublets()
        if should_plot:
            logger.info("Plotting clone graph")
            clone_graph.plot(output_dir / "graph", highlight_cell_ids, doublets2)

        logger.info(f"Removing {len(doublets2)} doublets from the graph (second round)")
        clone_graph.remove_nodes(doublets2)

        with open(output_dir / "doublets.txt", "w") as doublets_file:
            for clone in doublets + doublets2:
                assert clone.n == 1
                print(clone.cell_ids[0], file=doublets_file)

    if should_plot:
        logger.info("Plotting corrected clone graph")
        clone_graph.plot(output_dir / "graph_corrected", highlight_cell_ids)

    with open(output_dir / "components_corrected.txt", "w") as components_file:
        print(
            clone_graph.components_txt(highlight_cell_ids), file=components_file, end=""
        )

    clones = clone_graph.clones()
    with open(output_dir / "clones.txt", "w") as f:
        clone_graph.write_clones(f, clones)
    with open(output_dir / "clone_sequences.txt", "w") as f:
        clone_graph.write_clone_sequences(f, clones)
    logger.info(f"Detected {len(clones)} clones")
    clone_sizes = Counter(len(cells) for clone_id, cells, _ in clones)
    logger.info(
        "Clone size histogram\n size count\n%s",
        "\n".join(f"{k:5d} {clone_sizes[k]:5d}" for k in sorted(clone_sizes)),
    )
    unique_sizes = Counter(unique for _, _, unique in clones)
    logger.info(
        "Clone size histogram counting only unique cloneID combinations\n size count\n%s",
        "\n".join(f"{k:5d} {unique_sizes[k]:5d}" for k in sorted(unique_sizes)),
    )

    number_of_cells_in_clones = sum(k * v for k, v in clone_sizes.items())
    logger.debug("No. of cells in clones: %d", number_of_cells_in_clones)

    if should_write_loom:
        if len(transcriptome_inputs) > 1:
            logger.warning(
                "Writing a loom file only for the first transcriptome dataset"
            )
        cellranger = make_cellranger(transcriptome_inputs[0])
        write_loom(cells, cellranger, output_dir, clone_id_length=end - start)


def read_allowed_cellids(path):
    """
    Read a user-provided list of allowed cell IDs from a CSV

    Example:

    1   ACGTACGTACGTACGT_10x99

    or:

    1   ACGTACGTACGTACGT
    """
    allowed_ids = []
    filtered_df = pd.read_csv(Path(path), sep="\t", index_col=0, header=None)
    for cell_id in filtered_df.iloc[:, 0]:
        if cell_id.endswith("-1"):
            raise TrexError(
                "Cell ids in the list of allowed cell IDs must not end in '-1'"
            )
        allowed_ids.append(cell_id)
    logger.info(f"Restricting analysis to {len(allowed_ids)} allowed cells")
    return set(allowed_ids)


def read_excluded_clone_ids(path: Path) -> List[str]:
    """
    Read a user-provided list of cloneIDs to be ignored from a text file
    """
    excluded_clone_ids = pd.read_table(path, header=None)
    excluded_clone_ids = excluded_clone_ids[excluded_clone_ids.columns[0]].values
    logger.info(
        f"{len(excluded_clone_ids)} CloneIDs will be ignored during the analysis"
    )
    return set(excluded_clone_ids)


def is_similar(s: str, t: str, min_overlap: int, max_hamming: int) -> bool:
    if len(s) != len(t):
        raise IndexError("Sequences do not have the same length")

    matches = 0
    mismatches = 0
    for ch1, ch2 in zip(s, t):
        if ch1 == "-" or ch1 == "0" or ch2 == "-" or ch2 == "0":
            continue
        if ch1 == ch2:
            matches += 1
        else:
            mismatches += 1
    if matches + mismatches < min_overlap:
        return False
    if mismatches > max_hamming:
        return False
    return True
    # TODO allowed Hamming distance should be reduced relative to the overlap length


class SimilaritySet:
    def __init__(self, strings: set[str]):
        self._length = len(next(iter(strings)))
        if not all(len(s) == self._length for s in strings):
            raise ValueError("All strings must have the same length")
        self._exact = {s for s in strings if "0" not in s and "-" not in s}
        self._partial = [s for s in strings if "0" in s or "-" in s]

    def _is_similar(self, s, t):
        for ch1, ch2 in zip(s, t):
            if ch1 == "-" or ch1 == "0" or ch2 == "-" or ch2 == "0":
                continue
            if ch1 != ch2:
                return False
        return True

    def _occurs_partially(self, s: str) -> bool:
        for t in self._exact:
            if self._is_similar(s, t):
                return True
        return False

    def contains(self, s: str) -> bool:
        if s in self._exact:
            return True
        for t in self._partial:
            if is_similar(s, t, 0, 0):
                return True
        if ("0" in s or "-" in s) and self._occurs_partially(s):
            return True
        return False


def correct_clone_ids(
    molecules: List[Molecule], max_hamming: int, min_overlap: int = 20
) -> List[Molecule]:
    """
    Attempt to correct sequencing errors in the cloneID sequences of all molecules
    """
    # Obtain all cloneIDs (including those with '-' and '0')
    clone_ids = [m.clone_id for m in molecules]

    # Count the full-length cloneIDs
    counts = Counter(clone_ids)

    # Cluster them by Hamming distance
    clusters = cluster_sequences(
        list(set(clone_ids)),
        is_similar=lambda s, t: is_similar(s, t, min_overlap, max_hamming),
        k=7,
    )

    # Map non-singleton cloneIDs to a cluster representative
    clone_id_map = dict()
    for cluster in clusters:
        if len(cluster) > 1:
            # Pick most frequent cloneID as representative
            representative = max(cluster, key=lambda bc: (counts[bc], bc))
            for clone_id in cluster:
                clone_id_map[clone_id] = representative

    # Create a new list of molecules in which the cloneIDs have been replaced
    # by their representatives
    new_molecules = []
    for molecule in molecules:
        clone_id = clone_id_map.get(molecule.clone_id, molecule.clone_id)
        molecule = dataclasses.replace(molecule, clone_id=clone_id)
        new_molecules.append(molecule)
    return new_molecules


def correct_clone_ids_per_cell(
    molecules: List[Molecule], max_hamming: int, min_overlap: int = 20
) -> List[Molecule]:
    """
    Attempt to correct sequencing errors in the CloneID sequences of all
    molecules looking for similar sequences in the same cell.
    """
    # Count all cloneIDs (including those with '-' and '0')
    counts = Counter(m.clone_id for m in molecules)

    # Group molecules into cells
    cells: DefaultDict[str, List[Molecule]] = defaultdict(list)
    for molecule in molecules:
        cells[molecule.cell_id].append(molecule)

    # Iterate over cells
    cell_correction_map = defaultdict(dict)
    for cell_id, cell_molecules in cells.items():
        cell_clone_ids = list(set(m.clone_id for m in cell_molecules))

        if len(cell_clone_ids) < 2:
            continue
        clusters = cluster_sequences(
            cell_clone_ids,
            is_similar=lambda s, t: is_similar(s, t, min_overlap, max_hamming),
            k=0,
        )
        for cluster in clusters:
            if len(cluster) < 2:
                continue

            # Pick most frequent cloneID as representative
            longest = max(len(re.sub("[0-]", "", x)) for x in cluster)
            subcluster = [x for x in cluster if len(re.sub("[0-]", "", x)) == longest]
            representative = max(
                subcluster, key=lambda clone_id: (counts[clone_id], clone_id)
            )
            cell_correction_map[cell_id].update(
                {
                    clone_id: representative
                    for clone_id in cluster
                    if clone_id != representative
                }
            )

    def corrected_molecule(molecule):
        this_correction_map = cell_correction_map.get(molecule.cell_id, None)
        if this_correction_map is not None:
            new_clone_id = this_correction_map.get(molecule.clone_id, molecule.clone_id)
            molecule = dataclasses.replace(molecule, clone_id=new_clone_id)
        return molecule

    # Create a new list of molecules in which the cloneIDs have been replaced
    # by their representatives
    return [corrected_molecule(molecule) for molecule in molecules]


def filter_visium(
    cells: Iterable[Cell],
    molecules: Iterable[Molecule],
) -> List[Cell]:
    """
    Filter: cloneIDs that have only a count of one and are also only based on one read are  removed
    """
    new_cells = []
    del_cells = 0
    del_cloneids = 0
    for cell in cells:
        cell_id = cell.cell_id
        counts = cell.counts.copy()
        for clone_id, count in cell.counts.items():
            if count > 1:
                # This cloneID occurs more than once in this cell - keep it
                continue
            for molecule in molecules:
                if (
                    molecule.cell_id == cell_id
                    and molecule.clone_id == clone_id
                    and molecule.read_count == 1
                ):
                    # This cloneID has only a read count of 1 - remove it
                    del_cloneids += 1
                    del counts[clone_id]
        if counts:
            new_cells.append(Cell(cell_id=cell.cell_id, counts=counts))
        else:
            del_cells += 1

    logger.info(
        f"Found {del_cloneids} single-read cloneIDs and removed {del_cells} cells"
    )
    return new_cells


def filter_cells(
    cells: Iterable[Cell],
    molecules: Iterable[Molecule],
    keep_single_reads: bool = False,
) -> List[Cell]:
    """
    Filter cloneIDs according to two criteria:

    - CloneIDs that have only a count of one and can be found in another cell are most
      likely results of contamination and are removed,
    - If keep_single_reads is False, cloneIDs that have only a count of one and are also only based
      on one read are also removed
    """
    overall_counts: Dict[str, int] = Counter()
    for cell in cells:
        overall_counts.update(cell.counts)

    single_read_clone_ids = set()
    for molecule in molecules:
        if molecule.read_count == 1:
            single_read_clone_ids.add(molecule.clone_id)
    logger.info(f"Found {len(single_read_clone_ids)} single-read cloneIDs")

    # Filter out cloneIDs with a count of one that appear in another cell
    new_cells = []
    for cell in cells:
        counts = cell.counts.copy()
        for clone_id, count in cell.counts.items():
            if count > 1:
                # This cloneID occurs more than once in this cell - keep it
                continue
            if overall_counts[clone_id] > 1:
                # This cloneID occurs also in other cells - remove it
                del counts[clone_id]
            elif clone_id in single_read_clone_ids and not keep_single_reads:
                del counts[clone_id]
        if counts:
            new_cells.append(Cell(cell_id=cell.cell_id, counts=counts))
    return new_cells

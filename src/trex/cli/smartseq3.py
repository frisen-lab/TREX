"""
Run on Smart-seq3 data processed by zUMIs software.
"""

import sys
import logging
from pathlib import Path
from collections import Counter
from typing import List, Dict, Iterable, Optional

from .run10x import read_allowed_cellids, correct_clone_ids
from . import (
    CommandLineError,
    add_file_logging,
    make_output_dir,
    add_common_arguments,
)
from .. import __version__
from ..writers import write_count_matrix, write_cells, write_reads_or_molecules
from ..clone import CellGraph
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
        run_smartseq3(
            output_dir,
            genome_name=args.genome_name,
            allowed_cell_ids=allowed_cell_ids,
            chromosome=args.chromosome,
            start=args.start - 1 if args.start is not None else None,
            end=args.end,
            transcriptome_inputs=transcriptome_inputs,
            amplicon_inputs=amplicon_inputs,
            sample_names=sample_names,
            prefix=args.prefix,
            max_hamming=args.max_hamming,
            min_length=args.min_length,
            jaccard_threshold=args.jaccard_threshold,
            keep_single_reads=args.keep_single_reads,
            should_write_umi_matrix=args.umi_matrix,
            should_plot=args.plot,
            highlight_cell_ids=highlight_cell_ids,
        )
    except TrexError as e:
        raise CommandLineError("%s", e)


def add_arguments(parser):
    groups = add_common_arguments(parser, smartseq=True)

    groups.filter.add_argument(
        "--keep-single-reads",
        action="store_true",
        default=False,
        help="Keep cloneIDs supported by only a single read. Default: Discard them",
    )

    groups.output.add_argument(
        "--umi-matrix",
        default=False,
        action="store_true",
        help="Create a UMI count matrix 'umi_count_matrix.csv' with "
        "cells as columns and cloneIDs as rows",
    )


def run_smartseq3(
    output_dir: Path,
    *,
    transcriptome_inputs: List[Path],
    amplicon_inputs: List[Path],
    genome_name: Optional[str] = None,
    allowed_cell_ids: Optional[List[str]] = None,
    chromosome: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    sample_names: Optional[List[str]] = None,
    prefix: bool = False,
    max_hamming: int = 5,
    min_length: int = 20,
    jaccard_threshold: float = 0.7,
    keep_single_reads: bool = False,
    should_write_umi_matrix: bool = False,
    should_plot: bool = False,
    highlight_cell_ids: Optional[List[str]] = None,
):
    if sample_names is not None and len(sample_names) != len(set(sample_names)):
        raise TrexError("The sample names need to be unique")

    dataset_reader = DatasetReader(
        output_dir, genome_name, chromosome, start, end, prefix
    )
    reads = dataset_reader.read_all(
        transcriptome_inputs,
        amplicon_inputs,
        sample_names,
        allowed_cell_ids,
        require_umis=True,
        cell_id_tag="BC",
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

    cells = compute_cells(corrected_molecules, min_length)
    logger.info(f"Detected {len(cells)} cells")
    write_cells(output_dir / "cells.txt", cells)

    cells = filter_cells(cells, corrected_molecules, keep_single_reads)
    logger.info(f"{len(cells)} filtered cells remain")
    write_cells(output_dir / "cells_filtered.txt", cells)

    if should_write_umi_matrix:
        logger.info("Writing UMI matrix")
        write_count_matrix(output_dir / "umi_count_matrix.csv", cells)

    clone_graph = CellGraph(cells, jaccard_threshold=jaccard_threshold)

    with open(output_dir / "components.txt", "w") as components_file:
        print(
            clone_graph.components_txt(highlight_cell_ids), file=components_file, end=""
        )
    if should_plot:
        logger.info("Plotting clone graph")
        clone_graph.plot(output_dir / "graph", highlight_cell_ids)

    bridges = clone_graph.bridges()
    logger.info(f"Removing {len(bridges)} bridges from the graph")
    clone_graph.remove_edges(bridges)
    with open(output_dir / "components_corrected.txt", "w") as components_file:
        print(
            clone_graph.components_txt(highlight_cell_ids), file=components_file, end=""
        )

    if should_plot:
        logger.info("Plotting corrected clone graph")
        clone_graph.plot(output_dir / "graph_corrected", highlight_cell_ids)

    clones = clone_graph.clones()
    with open(output_dir / "clones.txt", "w") as f:
        clone_graph.write_clones(f, clones)
    with open(output_dir / "clone_sequences.txt", "w") as f:
        clone_graph.write_clone_sequences(f, clones)
    logger.info(f"Detected {len(clones)} clones")
    clone_sizes = Counter(len(cells) for clone_id, cells in clones)
    logger.info(
        "Clone size histogram\n size count\n%s",
        "\n".join(f"{k:5d} {clone_sizes[k]:5d}" for k in sorted(clone_sizes)),
    )
    number_of_cells_in_clones = sum(k * v for k, v in clone_sizes.items())
    logger.debug("No. of cells in clones: %d", number_of_cells_in_clones)
    assert len(cells) == number_of_cells_in_clones


def filter_cells(
    cells: Iterable[Cell],
    molecules: Iterable[Molecule],
    keep_single_reads: bool = False,
) -> List[Cell]:
    """
    Filter cloneIDs according to two criteria:

    - CloneIDs that have only a count of one and can be found in another cell are most
      likely results of contamination and are removed,
    - If keep_single_reads is False, cloneIDs that have only a count of one and are
      also only based on one read are also removed
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

"""
Run on single cell 10X Chromium or spatial Visium data processed by Cell / Space Ranger software
"""
import sys
import logging
from pathlib import Path
from collections import Counter
from typing import List, Dict, Iterable

import pandas as pd

from . import (
    setup_logging,
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
from ..clone import CloneGraph
from ..molecule import (Molecule, 
                        compute_molecules, 
                        correct_clone_ids, 
                        correct_clone_ids_per_cell, 
                        remove_odd_clone_ids)
from ..cell import Cell, compute_cells
from ..error import TrexError
from ..dataset import DatasetReader


__author__ = "leonie.von.berlin@ki.se"

logger = logging.getLogger(__name__)


def main(args):
    setup_logging(debug=args.debug)

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
        run_trex(
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
            min_bases_detected=args.min_bases_detected,
            per_cell_correction=args.per_cell,
            max_hamming=args.max_hamming,
            min_length=args.min_length,
            jaccard_threshold=args.jaccard_threshold,
            keep_single_reads=args.keep_single_reads,
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
        "--visium",
        default=False,
        action="store_true",
        help="Adjust filter settings for 10x Visium data: Filter out cloneIDs only based on "
        "one read, but keep those with only one UMI",
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
    genome_name: str,
    allowed_cell_ids: List[str],
    chromosome: str,
    start: int,
    end: int,
    transcriptome_inputs: List[Path],
    amplicon_inputs: List[Path],
    sample_names: List[str],
    prefix: bool,
    min_bases_detected: int,
    per_cell_correction: bool,
    max_hamming: int,
    min_length: int,
    jaccard_threshold: float,
    keep_single_reads: bool,
    should_write_umi_matrix: bool,
    should_run_visium: bool,
    should_plot: bool,
    highlight_cell_ids: List[str],
    should_write_loom: bool,
):
    if len(sample_names) != len(set(sample_names)):
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

    molecules = remove_odd_barcodes(molecules, min_bases_detected)

    write_reads_or_molecules(output_dir / "molecules.txt", molecules,
                             sort=False)

    if per_cell_correction:
        corrected_molecules = correct_barcodes_per_cell(molecules, max_hamming,
                                                        min_length)
    else:
        corrected_molecules = correct_clone_ids(molecules, max_hamming,
                                                min_length)
    
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
            clone_graph.components_txt(highlight_cell_ids), 
            file=components_file, end="",
        )
    if should_plot:
        logger.info("Plotting clone graph")
        clone_graph.plot(output_dir / "graph", highlight_cell_ids)

    bridges = clone_graph.bridges()
    logger.info(f"Removing {len(bridges)} bridges from the graph")
    clone_graph.remove_edges(bridges)
    with open(output_dir / "components_corrected.txt", "w") as components_file:
        print(
            clone_graph.components_txt(highlight_cell_ids), 
            file=components_file, end="",
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

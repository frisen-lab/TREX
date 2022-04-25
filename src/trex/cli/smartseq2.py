"""
Run on Smart-Seq2 data
"""
import sys
import logging
from pathlib import Path
from collections import Counter
from typing import List, Iterable

from .run10x import read_allowed_cellids, correct_clone_ids
from . import setup_logging, CommandLineError, add_file_logging, make_output_dir
from .. import __version__
from ..writers import write_count_matrix, write_cells, write_reads_or_molecules
from ..clone import CloneGraph
from ..molecule import Molecule
from ..error import TrexError
from ..cell import Cell, compute_cells
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
            f'Output directory "{output_dir}" already exists '
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
        run_smartseq2(
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
            readcount_threshold=args.readcount_threshold,
            should_write_read_matrix=args.read_matrix,
            should_plot=args.plot,
            highlight_cell_ids=highlight_cell_ids,
        )
    except TrexError as e:
        raise CommandLineError("%s", e)


def add_arguments(parser):
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument(
        "--debug",
        default=False,
        action="store_true",
        help="Print some extra debugging messages",
    )
    parser.add_argument(
        "--genome-name",
        metavar="NAME",
        help="Name of the genome as indicated in Cell Ranger count run with the flag --genome. "
        "Default: None",
        default=None,
    )
    parser.add_argument(
        "--chromosome",
        "--chr",
        help="Name of chromosome on which clone ID is located."
        " Default: Last chromosome in BAM file",
        default=None,
    )
    parser.add_argument(
        "--output",
        "-o",
        "--name",
        "-n",
        metavar="DIRECTORY",
        type=Path,
        help="name of the run and directory created by program. Default: %(default)s",
        default=Path("trex_run"),
    )
    parser.add_argument(
        "--delete", action="store_true", help="Delete output directory if it exists"
    )
    parser.add_argument(
        "--start",
        "-s",
        help="Position of first clone ID nucleotide. Default: Auto-detected",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--end",
        "-e",
        help="Position of last clone ID nucleotide. Default: Auto-detected",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--min-length",
        "-m",
        help="Minimum number of nucleotides a clone ID must have. Default: %(default)s",
        type=int,
        default=20,
    )
    parser.add_argument(
        "--max-hamming",
        help="Maximum hamming distance allowed for two clone IDs to be called similar. "
        "Default: %(default)s",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--jaccard-threshold",
        type=float,
        default=0,
        metavar="VALUE",
        help="If the Jaccard index between clone IDs of two cells is higher than VALUE, the cells "
        "are considered similar. Default: %(default)s",
    )
    parser.add_argument(
        "--amplicon",
        "-a",
        nargs="+",
        metavar="DIRECTORY",
        help="Path to united bam file for all cells or path to a folder with one bam file per cell,"
        "containing sequencing of the clone ID amplicon library. Provide these in "
        "same order as transcriptome datasets",
        default=None,
    )
    parser.add_argument(
        "--filter-cellids",
        "-f",
        metavar="CSV",
        type=Path,
        help="CSV file containing cell IDs to keep in the analysis."
        " This flag enables to remove cells e.g. doublets",
        default=None,
    )
    parser.add_argument(
        "--readcount-threshold",
        default=2,
        type=int,
        help="Minimum number of reads supporting a clone ID in order to keep it for downstream analysis. "
        "Default: %(default)s",
    )
    parser.add_argument(
        "--highlight", help="Highlight cell IDs listed in FILE in the clone graph"
    )
    parser.add_argument(
        "--samples",
        help="Sample names separated by comma, in the same order as bam files/bam file directories",
        default=None,
    )
    parser.add_argument(
        "--prefix",
        default=False,
        action="store_true",
        help="Add sample name as prefix to cell IDs (instead of as suffix)",
    )
    parser.add_argument(
        "--read-matrix",
        default=False,
        action="store_true",
        help="Creates a read count matrix with cells as columns and clone IDs as rows",
    )
    parser.add_argument(
        "--plot",
        dest="plot",
        default=False,
        action="store_true",
        help="Plot the clone graph",
    )
    parser.add_argument(
        "path",
        type=Path,
        nargs="+",
        metavar="DIRECTORY",
        help="Path to a united bam file for all cells or path to a folder with one bam file per cell.",
    )


def run_smartseq2(
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
    max_hamming: int,
    min_length: int,
    jaccard_threshold: float,
    readcount_threshold: int,
    should_write_read_matrix: bool,
    should_plot: bool,
    highlight_cell_ids: List[str],
):
    if len(sample_names) != len(set(sample_names)):
        raise TrexError("The sample names need to be unique")

    dataset_reader = DatasetReader(
        output_dir, genome_name, chromosome, start, end, prefix
    )
    reads = dataset_reader.read_all(
        transcriptome_inputs,
        amplicon_inputs,
        sample_names,
        allowed_cell_ids,
        require_umis=False,
        cell_id_tag="BC",
    )
    if not reads:
        raise TrexError("No reads left after --filter-cellids filtering")

    clone_ids = [
        r.clone_id for r in reads if "-" not in r.clone_id and "0" not in r.clone_id
    ]
    logger.info(
        f"Read {len(reads)} reads containing (parts of) the clone ID "
        f"({len(clone_ids)} full clone IDs, {len(set(clone_ids))} unique)"
    )

    write_reads_or_molecules(output_dir / "reads.txt", reads, require_umis=False)

    # We do not have multiple reads per molecule, so we treat each read as one molecule
    molecules = reads

    corrected_molecules = correct_clone_ids(molecules, max_hamming, min_length)
    clone_ids = [
        m.clone_id
        for m in corrected_molecules
        if "-" not in m.clone_id and "0" not in m.clone_id
    ]
    logger.info(
        f"After clone ID correction, {len(set(clone_ids))} unique clone IDs remain"
    )

    write_reads_or_molecules(
        output_dir / "molecules_corrected.txt",
        corrected_molecules,
        require_umis=False,
        sort=False,
    )

    cells = compute_cells(corrected_molecules, min_length)
    logger.info(f"Detected {len(cells)} cells")
    write_cells(output_dir / "cells.txt", cells)

    cells = filter_smartseq(cells, corrected_molecules, readcount_threshold)
    logger.info(f"{len(cells)} filtered cells remain")
    write_cells(output_dir / "cells_filtered.txt", cells)

    if should_write_read_matrix:
        logger.info("Writing read matrix")
        write_count_matrix(output_dir / "read_count_matrix.csv", cells)

    clone_graph = CloneGraph(cells, jaccard_threshold=jaccard_threshold)

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


def filter_smartseq(
    cells: Iterable[Cell], molecules: Iterable[Molecule], readcount_threshold: int
) -> List[Cell]:
    """
    Removes clone IDs that are supported by less reads than the given readcount_threshold
    """
    new_cells = []
    for cell in cells:
        counts = cell.counts.copy()
        for clone_id, count in cell.counts.items():
            if count < readcount_threshold:
                del counts[clone_id]
        if counts:
            new_cells.append(Cell(cell_id=cell.cell_id, counts=counts))
    return new_cells

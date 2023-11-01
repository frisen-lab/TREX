import logging
from pathlib import Path
import shutil
from types import SimpleNamespace

logger = logging.getLogger(__name__)


class CommandLineError(Exception):
    pass


class NiceFormatter(logging.Formatter):
    """
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).

    Based on http://stackoverflow.com/a/9218261/715090 .
    """

    def format(self, record):
        if record.levelno != logging.INFO:
            record.msg = "{}: {}".format(record.levelname, record.msg)
        return super().format(record)


def setup_logging(debug: bool) -> None:
    """
    Set up logging. If debug is True, then DEBUG level messages are printed.
    """
    handler = logging.StreamHandler()
    handler.setFormatter(NiceFormatter())

    root = logging.getLogger()
    root.addHandler(handler)
    root.setLevel(logging.DEBUG if debug else logging.INFO)


def add_file_logging(path: Path) -> None:
    file_handler = logging.FileHandler(path)
    root = logging.getLogger()
    root.addHandler(file_handler)


def make_output_dir(path, delete_if_exists):
    try:
        path.mkdir()
    except FileExistsError:
        if delete_if_exists:
            logger.debug(f'Re-creating folder "{path}"')
            shutil.rmtree(path)
            path.mkdir()
        else:
            raise


def add_common_arguments(parser, smartseq: bool):
    """Add arguments to an ArgumentParser common to both run10x and smartseq2/3"""

    input_group = parser.add_argument_group("Input")

    input_group.add_argument(
        "--genome-name",
        metavar="NAME",
        help="Name of the genome as indicated in 'cellranger count' run with the flag --genome. "
        "Default: Auto-detected",
        default=None,
    )
    input_group.add_argument(
        "--chromosome",
        "--chr",
        help="Name of chromosome on which cloneID is located. "
        "Default: Last chromosome in BAM file",
        default=None,
    )
    input_group.add_argument(
        "--start",
        "-s",
        help="Position of first cloneID nucleotide (1-based). Default: Auto-detected",
        type=int,
        metavar="INT",
        default=None,
    )
    input_group.add_argument(
        "--end",
        "-e",
        help="Position of last cloneID nucleotide (1-based). Default: Auto-detected",
        type=int,
        metavar="INT",
        default=None,
    )
    if smartseq:
        help = (
            "Path to united BAM file for all cells or path to a folder with one BAM file "
            "per cell, containing sequencing of the cloneID amplicon library."
        )
    else:
        help = (
            "Path to Cell Ranger result directory (a subdirectory 'outs' must exist) "
            "containing sequencing of the cloneID amplicon library."
        )
    input_group.add_argument(
        "--amplicon",
        "-a",
        nargs="+",
        metavar="DIRECTORY",
        help=help + " Provide these in the same order as transcriptome datasets",
        default=None,
    )
    if smartseq:
        help = "BAM files/BAM file directories"
    else:
        help = "Cell Ranger directories"
    input_group.add_argument(
        "--samples",
        help="Sample names separated by comma, in the same order as " + help,
        default=None,
    )
    input_group.add_argument(
        "--prefix",
        default=False,
        action="store_true",
        help="Add sample name as prefix to cell IDs. Default: Add as suffix",
    )

    filter_group = parser.add_argument_group("Filter settings")

    filter_group.add_argument(
        "--min-length",
        "-m",
        help="Minimum number of nucleotides a cloneID must have. Default: %(default)s",
        type=int,
        metavar="INT",
        default=20,
    )
    filter_group.add_argument(
        "--max-hamming",
        help="Maximum hamming distance allowed for two cloneIDs to be called similar. "
        "Default: %(default)s",
        type=int,
        metavar="INT",
        default=5,
    )
    filter_group.add_argument(
        "--jaccard-threshold",
        type=float,
        default=0.7,
        metavar="VALUE",
        help="If the Jaccard index between cloneIDs of two cells is higher than VALUE, they "
        "are considered similar. Default: %(default)s",
    )
    filter_group.add_argument(
        "--filter-cellids",
        "-f",
        metavar="TSV",
        type=Path,
        help="TSV file containing cell IDs to keep in the analysis. "
        "This flag enables to remove cells e.g. doublets. "
        "Expected format: see documentation",
        default=None,
    )
    filter_group.add_argument(
        "--filter-cloneids",
        type=Path,
        help="Text file with cloneIDs to be ignored during the analysis. Format: One cloneID per line. "
        "Use this to remove e.g. overrepresented cloneIDs or misalignments.",
    )

    output_group = parser.add_argument_group("Output directory")

    output_group.add_argument(
        "--output",
        "-o",
        "--name",
        "-n",
        metavar="DIRECTORY",
        type=Path,
        help="Name of the run directory to be created by the program. Default: %(default)s",
        default=Path("trex_run"),
    )
    output_group.add_argument(
        "--delete",
        action="store_true",
        help="Delete the run directory if it already exists",
    )

    optional_group = parser.add_argument_group(
        "Optional output files",
        description="Use these options to enable creation "
        "of additional files in the output directory",
    )
    optional_group.add_argument(
        "--plot",
        dest="plot",
        default=False,
        action="store_true",
        help="Plot the clone graph. This requires GraphViz to be installed.",
    )
    optional_group.add_argument(
        "--highlight",
        metavar="FILE",
        help="Highlight cell IDs listed in FILE "
        "(text file with one cell ID per line) in the clone graph",
    )

    if smartseq:
        help = (
            "Path to a united BAM file for all cells or path to a folder "
            "with one BAM file per cell."
        )
    else:
        help = (
            "Path to the input Cell Ranger directories. "
            "There must be an 'outs' subdirectory in each of these directories."
        )
    parser.add_argument(
        "path",
        type=Path,
        nargs="+",
        metavar="DIRECTORY",
        help=help,
    )

    return SimpleNamespace(input=input_group, filter=filter_group, output=output_group)

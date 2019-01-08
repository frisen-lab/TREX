"""
Extract and filter random barcodes from single-cell sequencing data
"""
import sys
import argparse
from collections import Counter, defaultdict, OrderedDict
import operator
import warnings
import logging
from pathlib import Path
from typing import Set, List, Dict, NamedTuple, Iterable, Callable

import numpy as np
import pysam
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', 'Conversion of the second argument of issubdtype')
    import loompy

__version__ = '0.1'
__author__ = 'leonie.von.berlin@stud.ki.se'

logger = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('--genome-name', metavar='NAME',
        help='Name of the genome as indicated in cell ranger count run with the flag --genome. '
             'Default: Auto-detected',
        default=None)
    parser.add_argument('--chromosome', '--chr',
        help='Barcode chromosome name. Default: Last chromosome in the BAM file',
        default=None)
    parser.add_argument('--output', '-o', '--name', '-n', metavar='DIRECTORY', type=Path,
        help='name of the run and directory created by program. Default: %(default)s',
        default=Path('lineage_run'))
    parser.add_argument('--start', '-s',
        help='Position of first barcode base. Default: Auto-detected',
        type=int, default=None)
    parser.add_argument('--end', '-e',
        help='Position of last barcode base. Default: Auto-detected',
        type=int, default=None)
    parser.add_argument('--min-length', '-m',
        help='Minimum number of bases a barcode must have. Default: %(default)s',
        type=int, default=20)
    parser.add_argument('--max-hamming',
        help='Hamming distance allowed for two barcodes to be called similar. '
            'Default: %(default)s',
        type=int, default=5)
    parser.add_argument('--keep-single-reads', action='store_true', default=False,
        help='Keep barcodes supported by only a single read. Default: Discard them')
    parser.add_argument('-l', '--loom',
        help='If given, create loom-file from cell ranger and barcode data. '
            'File will have the same name as the run',
        action='store_true')
    parser.add_argument('path', metavar='DIRECTORY', type=Path,
        help='Path to cell ranger "outs" directory')
    return parser.parse_args()


def hamming_distance(s, t):
    """Return Hamming distance between s and t."""
    # This explicit for loop is slightly faster than
    # using the comprehension sum(1 for c, d in zip(s, t) if c != d)
    n = 0
    for c, d in zip(s, t):
        if c != d:
            n += 1
    return n


class Graph:
    """Graph that can find connected components"""

    def __init__(self, nodes):
        # values are lists of adjacent nodes
        self._nodes = OrderedDict()
        for node in nodes:
            self._nodes[node] = []

    def add_edge(self, node1, node2):
        self._nodes[node1].append(node2)
        self._nodes[node2].append(node1)

    def connected_components(self):
        """Return a list of connected components."""
        visited = set()
        components = []
        for node, neighbors in self._nodes.items():
            if node in visited:
                continue
            # Start a new component
            to_visit = [node]
            component = []
            while to_visit:
                n = to_visit.pop()
                if n in visited:
                    continue
                visited.add(n)
                component.append(n)
                for neighbor in self._nodes[n]:
                    if neighbor not in visited:
                        to_visit.append(neighbor)
            components.append(component)
        return components


def kmers(s: str, k: int):
    """
    Yield all overlapping k-kmers of s.

    >>> list(kmers('hello', 3))
    ['hel', 'ell', 'llo']
    """
    for i in range(len(s) - k + 1):
        yield s[i:i+k]


def cluster_sequences(
        sequences: List[str], is_similar: Callable[[str, str], bool], k: int=6) -> List[List[str]]:
    """
    Cluster sequences by Hamming distance.

    If k > 0, a k-mer filter may be enabled for speedups. With the filter, a single barcode is not
    compared to all others, but only to those others with which it shares a k-mer.

    Speedups happen only at k >= 5, so for values lower than that, the filter is disabled.

    The filter is a heuristic, so results may differ when it is enabled. This is relevant
    starting at about k=8.

    Each element of the returned list is a cluster.
    """
    graph = Graph(sequences)
    if k >= 5:
        shared_kmers = defaultdict(set)
        for bc in sequences:
            for kmer in kmers(bc, k):
                shared_kmers[kmer].add(bc)

        for bc in sequences:
            others: Set[str] = set()
            for kmer in kmers(bc, k):
                others.update(shared_kmers[kmer])
            for other in others:
                if is_similar(bc, other):
                    graph.add_edge(bc, other)
    else:
        for i, x in enumerate(sequences):
            for j in range(i+1, len(sequences)):
                y = sequences[j]
                assert len(x) == len(y)
                if is_similar(x, y):
                    graph.add_edge(x, y)
    return graph.connected_components()


def read_cellids(path: Path) -> Set[str]:
    """
    Read barcodes.tsv, which contains a list of corrected and approved cellIDs like this:

    AAACCTGAGCGACGTA-1
    AAACCTGCATACTCTT-1
    """
    with open(path) as f:
        ids = []
        for line in f:
            line = line.strip('\n')
            ids.append(line)
    return set(ids)


class Read(NamedTuple):
    cell_id: str
    umi: str
    lineage_id: str


class Molecule(NamedTuple):
    cell_id: str
    umi: str
    lineage_id: str
    read_count: int


class Cell(NamedTuple):
    cell_id: str
    lineage_id_counts: Dict[str, int]


def detect_lineage_id_location(alignment_file, reference_name):
    """
    Detect where the lineage id is located on the reference by inspecting the alignments.

    Return (lineage_id_start, lineage_id_end)
    """
    # Look for reference positions at which reads are soft-clipped at their 3' end
    starts = Counter()
    reference_length = alignment_file.get_reference_length(reference_name)
    for alignment in alignment_file.fetch(reference_name):
        clip_right = alignment.query_length - alignment.query_alignment_end
        if clip_right >= 5:
            starts[alignment.reference_end] += 1

    for lineage_id_start, freq in starts.most_common(5):
        # Soft-clipping at the 5' end cannot be used to find the lineage id end when
        # the lineage id region is too far at the 3' end of the contig. Instead,
        # look at pileups and check base frequencies (over the lineage id, bases should
        # be roughly uniformly distributed).
        if lineage_id_start >= reference_length:
            # The most common reference position that is soft clipped is often the 3' end
            # of the contig. Skip that.
            continue
        lineage_id_end = lineage_id_start
        for column in alignment_file.pileup(reference_name, start=lineage_id_start):
            if column.reference_pos < lineage_id_start:
                # See pileup() documentation
                continue
            bases = [p.alignment.query_sequence[p.query_position] for p in column.pileups if p.query_position is not None]
            counter = Counter(bases)
            # Check whether one base dominates
            if counter.most_common()[0][1] / len(bases) > 0.95:
                # We appear to have found the end of the lineage id
                lineage_id_end = column.reference_pos
                break
        if lineage_id_end - lineage_id_start >= 5:
            # Good enough
            return (lineage_id_start, lineage_id_end)
    raise ValueError(f'Could not detect lineage id location on chromosome {reference_name}')


def read_bam(bam_path: Path, output_dir: Path, cell_ids, chr_name, lineage_id_start=None, lineage_id_end=None):
    """
    bam_path -- path to input BAM file
    output_bam_path -- path to an output BAM file. All reads on the chromosome that have the
        required tags are written to this file
    """
    with pysam.AlignmentFile(bam_path) as alignment_file:
        if chr_name is None:
            chr_name = alignment_file.references[-1]

        if lineage_id_start is None or lineage_id_end is None:
            if lineage_id_start is not None or lineage_id_end is not None:
                raise ValueError('Either both or none of lineage id start and end must be provided')
            lineage_id_start, lineage_id_end = detect_lineage_id_location(alignment_file, chr_name)
        logger.info(f'Reading lineage ids from {chr_name}:{lineage_id_start + 1}-{lineage_id_end}')
        if lineage_id_end - lineage_id_start < 10:
            raise ValueError('Auto-detected lineage id too short, something is wrong')
        output_bam_path = output_dir / (chr_name + '_entries.bam')
        with pysam.AlignmentFile(output_bam_path, 'wb', template=alignment_file) as out_bam:
            # Fetches those reads aligning to the artifical, lineage-id-containing chromosome
            reads = []
            unknown_ids = no_cell_id = no_umi = 0
            for read in alignment_file.fetch(chr_name):
                # Skip reads without cellID or UMI
                if not read.has_tag('CB'):
                    no_cell_id += 1
                if not read.has_tag('UB'):
                    no_umi += 1
                if not read.has_tag('CB') or not read.has_tag('UB'):
                    continue
                # Filters out reads that have not approved cellIDs
                cell_id = read.get_tag('CB')
                if cell_id not in cell_ids:
                    unknown_ids += 1
                    continue

                query_align_end = read.query_alignment_end
                query_align_start = read.query_alignment_start

                # Extract lineage id
                lineage_id = ['-'] * (lineage_id_end - lineage_id_start)
                bases = 0
                for query_pos, ref_pos in read.get_aligned_pairs():
                    # Replace soft-clipping with an ungapped alignment extending into the
                    # soft-clipped region, assuming the clipping occurred because the lineage id
                    # region was encountered.
                    if ref_pos is None:
                        # Soft clip or insertion
                        if query_align_end <= query_pos:
                            # We are in a soft-clipped region at the 3' end of the read
                            ref_pos = read.reference_end + (query_pos - query_align_end)
                        elif query_align_start > query_pos:
                            # We are in a soft-clipped region at the 5' end of the read
                            ref_pos = read.reference_start - (query_align_start - query_pos)
                        # ref_pos remains None if this is an insertion

                    if ref_pos is not None and lineage_id_start <= ref_pos < lineage_id_end:
                        if query_pos is None:
                            # Deletion or intron skip
                            query_base = '0'
                        else:
                            # Match or mismatch
                            query_base = read.query_sequence[query_pos]
                            bases += 1
                        lineage_id[ref_pos - lineage_id_start] = query_base
                if bases == 0:
                    # Skip if this read does not cover the lineage id
                    continue
                lineage_id = ''.join(lineage_id)
                reads.append(Read(cell_id=cell_id, umi=read.get_tag('UB'), lineage_id=lineage_id))

                # Write the passing alignments to a separate file
                out_bam.write(read)

            logger.info(f'Skipped {unknown_ids} reads with unrecognized cell ids '
                        f'(and {no_umi+no_cell_id} without UMI or cell id)')
    sorted_reads = sorted(reads, key=lambda read: (read.umi, read.cell_id, read.lineage_id))
    assert len(sorted_reads) == 0 or len(sorted_reads[0].lineage_id) == lineage_id_end - lineage_id_start
    return sorted_reads


def compute_consensus(sequences):
    """
    Compute a consensus for a set of sequences.

    All sequences must have the same length.
    """
    if len(sequences) == 1:
        return sequences[0]
    assert sequences

    # TODO
    # Ensure that the sequences are actually somewhat similar

    letters = np.array(['A', 'C', 'G', 'T', '-', '0'])
    length = len(sequences[0])
    consens_np = np.zeros([length, 6], dtype='float16')
    for sequence in sequences:
        align = np.zeros([length, 6], dtype='float16')
        for (i, ch) in enumerate(sequence):
            # turns each base into a number and position in numpy array
            if ch == 'A':
                align[i, 0] = 1
            elif ch == 'C':
                align[i, 1] = 1
            elif ch == 'G':
                align[i, 2] = 1
            elif ch == 'T':
                align[i, 3] = 1
            elif ch == '-':
                align[i, 4] = 0.1
            elif ch == '0':
                align[i, 5] = 0.1
        consens_np += align
    # calculate base with maximum count for each position
    bin_consens = np.argmax(align, axis=1)
    # convert maximum counts into consensus sequence
    return ''.join(letters[bin_consens])


def compute_molecules(sorted_reads):
    """
    - Forms groups of reads with identical CellIDs and UMIs => belong to one molecule

    - forms consensus sequence of all lineage ids of one group,
    """
    groups = defaultdict(list)
    for read in sorted_reads:
        groups[(read.umi, read.cell_id)].append(read.lineage_id)

    molecules = []
    for (umi, cell_id), lineage_ids in groups.items():
        lineage_id_consensus = compute_consensus(lineage_ids)
        molecules.append(
            Molecule(cell_id=cell_id, umi=umi, lineage_id=lineage_id_consensus,
                read_count=len(lineage_ids)))

    sorted_molecules = sorted(molecules, key=lambda mol: (mol.cell_id, mol.lineage_id, mol.umi))

    return sorted_molecules


def correct_lineage_ids(
        molecules: List[Molecule], max_hamming: int, min_overlap: int = 20) -> List[Molecule]:
    """
    Attempt to correct sequencing errors in the lineage id sequences of all molecules
    """
    # Obtain all lineage ids (including those with '-' and '0')
    lineage_ids = [m.lineage_id for m in molecules]

    # Count the full-length barcodes
    lineage_id_counts = Counter(lineage_ids)

    # Cluster them by Hamming distance
    def is_similar(s, t):
        # m = max_hamming
        if '-' in s or '-' in t:
            # Remove suffix and/or prefix where sequences do not overlap
            s = s.lstrip('-')
            t = t[-len(s):]
            s = s.rstrip('-')
            if len(s) < min_overlap:
                return False
            t = t[:len(s)]
            # TODO allowed Hamming distance should be reduced relative to the overlap length
            # m = max_hamming * len(s) / len(original_length_of_s)
        return hamming_distance(s, t) <= max_hamming

    clusters = cluster_sequences(list(set(lineage_ids)), is_similar=is_similar, k=6)

    # Map non-singleton barcodes to a cluster representative
    lineage_id_map = dict()
    for cluster in clusters:
        if len(cluster) > 1:
            # Pick most frequent lineage id as representative
            representative = max(cluster, key=lambda bc: (lineage_id_counts[bc], bc))
            for lineage_id in cluster:
                lineage_id_map[lineage_id] = representative

    # Create a new list of molecules in which the lineage ids have been replaced
    # by their representatives
    new_molecules = []
    for molecule in molecules:
        lineage_id = lineage_id_map.get(molecule.lineage_id, molecule.lineage_id)
        new_molecules.append(molecule._replace(lineage_id=lineage_id))
    return new_molecules


def compute_cells(sorted_molecules: List[Molecule], minimum_lineage_id_length: int) -> List[Cell]:
    """
    Group molecules into cells.
    """
    # 1. Forms groups of molecules (with set lineage id minimum length) that have identical cellIDs
    #    => belong to one cell,
    # 2. counts number of appearances of each lineage id in each group,

    cell_id_groups = defaultdict(list)
    for molecule in sorted_molecules:
        lineage_id = molecule.lineage_id
        pure_li = lineage_id.strip('-')
        # TODO may not work as intended (strip only removes prefixes and suffixes)
        pure_bc0 = lineage_id.strip('0')
        if len(pure_li) >= minimum_lineage_id_length and len(pure_bc0) >= minimum_lineage_id_length:
            cell_id_groups[molecule.cell_id].append(molecule)

    cells = []
    for cell_id, molecules in cell_id_groups.items():
        lineage_ids = [molecule.lineage_id for molecule in molecules]
        lineage_id_counts = OrderedDict(sorted(Counter(lineage_ids).most_common(),
            key=lambda x: x[0].count('-')))
        cells.append(Cell(cell_id=cell_id, lineage_id_counts=lineage_id_counts))
    return cells


def filter_cells(
        cells: Iterable[Cell], molecules: Iterable[Molecule],
        keep_single_reads: bool=False) -> List[Cell]:
    """
    Filter lineage ids according to two criteria:

    - Lineage ids that have only a count of one and can be found in another cell are most
      likely results of contamination and are removed,
    - If keep_single_reads is False, lineage ids that have only a count of one and are also only based
      on one read are also removed
    """
    overall_lineage_id_counts: Dict[str, int] = Counter()
    for cell in cells:
        overall_lineage_id_counts.update(cell.lineage_id_counts)

    single_read_lineage_ids = set()
    for molecule in molecules:
        if molecule.read_count == 1:
            single_read_lineage_ids.add(molecule.lineage_id)
    # or:
    # single_read_barcodes = {m.lineage_id for m in molecules if m.read_count == 1}

    # filters out lineage ids with a count of one that appear in another cell
    new_cells = []
    for cell in cells:
        lineage_id_counts = cell.lineage_id_counts.copy()
        for lineage_id, count in cell.lineage_id_counts.items():
            if count > 1:
                # This lineage id occurs more than once in this cell - keep it
                continue
            if overall_lineage_id_counts[lineage_id] > 1:
                # This lineage id occurs also in other cells - remove it
                del lineage_id_counts[lineage_id]
            elif lineage_id in single_read_lineage_ids and not keep_single_reads:
                del lineage_id_counts[lineage_id]
        new_cells.append(Cell(cell_id=cell.cell_id, lineage_id_counts=lineage_id_counts))
    return new_cells


def compute_lineages(cells: List[Cell]) -> Dict[str, List[Cell]]:
    lineage_id_counts: Dict[str, int] = Counter()
    for cell in cells:
        lineage_id_counts.update(cell.lineage_id_counts)
    all_lineage_ids = list(lineage_id_counts)

    # Create a graph of barcodes; add an edge for barcodes occuring in the same cell
    graph = Graph(all_lineage_ids)
    for cell in cells:
        barcodes = list(cell.lineage_id_counts)
        for other in barcodes[1:]:
            graph.add_edge(barcodes[0], other)
    clusters = graph.connected_components()

    cluster_sizes = Counter(len(c) for c in clusters)
    # logger.info(f'Lineage id cluster size histogram: {cluster_sizes}')

    # TODO copied from correct_lineage_ids
    # Map non-singleton barcodes to a cluster representative
    lineage_id_map = dict()
#    lineage_id_clusters = defaultdict(list)
    for cluster in clusters:
        # Pick most frequent lineage id as representative
        representative = max(cluster, key=lambda li: lineage_id_counts.get(li, 0))
        for lineage_id in cluster:
            lineage_id_map[lineage_id] = representative
        # Needs to be immutable as it will be shared by Lineage instances
#        lineage_id_clusters[representative] = tuple(cluster)

    # Group cells by the representative lineage id
    cell_groups = defaultdict(list)
    for cell in cells:
        if not cell.lineage_id_counts:
            continue
        # Since a cell can only belong to one lineage, only the first lineage id is used
        first_lineage_id = next(iter(cell.lineage_id_counts))
        representative = lineage_id_map[first_lineage_id]
        cell_groups[representative].append(cell)

    # TODO should we return List[Lineage]?
    return cell_groups


def write_loom(cells: List[Cell], cellranger_dir, output_dir, lineage_id_length, top_n=6):
    """
    Create a loom file from a Cell Ranger sample directory and augment it with information about
    the most abundant lineage id and their counts.
    """
    # For each cell, collect the most abundant lineage ids and their counts
    # Maps cell_id to a list of (lineage_id, count) pairs that represent the most abundant lineage ids.
    most_abundant = dict()
    for cell in cells:
        if not cell.lineage_id_counts:
            continue
        counts = sorted(cell.lineage_id_counts.items(), key=operator.itemgetter(1))
        counts.reverse()
        counts = counts[:top_n]
        most_abundant[cell.cell_id] = counts

    sample_dir = cellranger_dir.parent
    loompy.create_from_cellranger(sample_dir, outdir=output_dir)
    # create_from_cellranger() does not tell us the name of the created file,
    # so we need to re-derive it from the sample name.
    sample_name = sample_dir.name
    loom_path = output_dir / (sample_name + '.loom')

    with loompy.connect(loom_path) as ds:
        # Cell ids in the loom file are prefixed by the sample name and a ':'. Remove that prefix.
        loom_cell_ids = [cell_id[len(sample_name)+1:] for cell_id in ds.ca.CellID]

        # Transform lineage ids and count data
        # brings lineage id data into correct format for loom file.
        # Array must have same shape as all_cellIDs
        lineage_id_lists = [[] for _ in range(top_n)]
        count_lists = [[] for _ in range(top_n)]
        for cell_id in loom_cell_ids:
            lineage_id_counts = most_abundant.get(cell_id, [])
            # Fill up to a constant length
            while len(lineage_id_counts) < top_n:
                lineage_id_counts.append(('-', 0))

            for i, (lineage_id, count) in enumerate(lineage_id_counts):
                lineage_id_lists[i].append(lineage_id)
                count_lists[i].append(count)

        # Add lineage id and count information to loom file
        for i in range(top_n):
            ds.ca[f'linBarcode_{i+1}'] = np.array(lineage_id_lists[i], dtype='S%r' % lineage_id_length)
            ds.ca[f'linBarcode_count_{i+1}'] = np.array(count_lists[i], dtype=int)


def write_cells(path: Path, cells: List[Cell]) -> None:
    """Write cells to a tab-separated file"""
    with open(path, 'w') as f:
        print(
            '#Each output line corresponds to one cell and has the following style: '
            'CellID\t:\tBarcode1\tCount1\tBarcode2\tCount2...\n'
            '# dash (-) = barcode base outside of read, '
            '0 = deletion in barcode sequence (position unknown)', file=f)
        for cell in cells:
            row = [cell.cell_id, ':']
            sorted_lineage_ids = sorted(cell.lineage_id_counts, key=lambda x: cell.lineage_id_counts[x], reverse=True)
            if not sorted_lineage_ids:
                continue
            for lineage_id in sorted_lineage_ids:
                row.extend([lineage_id, cell.lineage_id_counts[lineage_id]])
            print(*row, sep='\t', file=f)


class NiceFormatter(logging.Formatter):
    """
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).

    Based on http://stackoverflow.com/a/9218261/715090 .
    """
    def format(self, record):
        if record.levelno != logging.INFO:
            record.msg = '{}: {}'.format(record.levelname, record.msg)
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
    file_handler.setFormatter(NiceFormatter())
    root = logging.getLogger()
    root.addHandler(file_handler)


def main():
    args = parse_arguments()
    input_dir = args.path
    output_dir = args.output
    setup_logging(debug=False)

    # PART I + II: Lineage id extraction and reads construction

    # 1. Extracts reads aligning to lineage id chromosome,
    # 2. extracts lineage ids, UMIs and cellIDs from reads,
    # 3. outputs UMI-sorted reads with barcodes

    try:
        output_dir.mkdir()
    except FileExistsError:
        logger.error(f'Output directory {output_dir!r} already exists '
            '(use -o to specify a different one)')
        sys.exit(1)

    add_file_logging(output_dir / 'log.txt')
    matrices_path = input_dir / 'filtered_gene_bc_matrices'
    if not matrices_path.exists():
        logger.error("Directory 'filtered_gene_bc_matrices/' must exist in the given path")
        sys.exit(1)
    if args.genome_name is None:
        genomes = [p for p in matrices_path.iterdir() if p.is_dir()]
        if not genomes:
            logger.error(f"No subfolders found in the 'outs/filtered_gene_bc_matrices/' folder")
            sys.exit(1)
        if len(genomes) > 1:
            logger.error('Exactly one genome folder expected in the '
                "'outs/filtered_gene_bc_matrices/' folder, but found:")
            for g in genomes:
                logger.error(f'  {g!r}')
            sys.exit(1)
        genome_dir = genomes[0]
    else:
        genome_dir = matrices_path / args.genome_name
    cell_ids = read_cellids(genome_dir / 'barcodes.tsv')
    logger.info(f'Found {len(cell_ids)} cell ids in the barcodes.tsv file')

    sorted_reads = read_bam(
        input_dir / 'possorted_genome_bam.bam', output_dir,
        cell_ids, args.chromosome, args.start - 1 if args.start is not None else None, args.end)

    lineage_ids = [
        r.lineage_id for r in sorted_reads if '-' not in r.lineage_id and '0' not in r.lineage_id]
    logger.info(f'Read {len(sorted_reads)} reads containing (parts of) the barcode '
        f'({len(lineage_ids)} full barcodes, {len(set(lineage_ids))} unique)')
    with open(output_dir / 'reads.txt', 'w') as reads_file:
        print(
            '#Each output line corresponds to one read and has the following style: '
            'CellID\tUMI\tBarcode\n'
            '# dash (-) = barcode base outside of read, '
            '0 = deletion in barcode sequence (position unknown)', file=reads_file)
        for read in sorted_reads:
            print(read.cell_id, read.umi, read.lineage_id, sep='\t', file=reads_file)

    # Part III: Molecule construction

    # 1. Forms groups of reads with identical CellIDs and UMIs => belong to one molecule,
    # 2. forms consensus sequence of all lineage ids of one group,
    # 3. outputs molecules and corresponding CellIDs/UMIs

    molecules = compute_molecules(sorted_reads)
    lineage_ids = [
        m.lineage_id for m in molecules if '-' not in m.lineage_id and '0' not in m.lineage_id]
    logger.info(f'Detected {len(molecules)} molecules ({len(lineage_ids)} full lineage ids, '
        f'{len(set(lineage_ids))} unique)')
    with open(output_dir / 'molecules.txt', 'w') as molecules_file:
        print(
            '#Each output line corresponds to one molecule and has the following style: '
            'CellID\tUMI\tBarcode\n'
            '# dash (-) = barcode base outside of read, '
            '0 = deletion in barcode sequence (position unknown)', file=molecules_file)
        for molecule in molecules:
            print(molecule.cell_id, molecule.umi, molecule.lineage_id, sep='\t', file=molecules_file)

    # Part IV: Cell construction

    corrected_molecules = correct_lineage_ids(molecules, args.max_hamming, args.min_length)
    lineage_ids = [m.lineage_id for m in corrected_molecules
        if '-' not in m.lineage_id and '0' not in m.lineage_id]
    logger.info(f'After lineage id correction, {len(set(lineage_ids))} unique lineage ids remain')

    with open(output_dir / 'molecules_corrected.txt', 'w') as molecules_file:
        print('cell_id', 'umi', 'lineage_id', sep='\t', file=molecules_file)
        for molecule in corrected_molecules:
            print(molecule.cell_id, molecule.umi, molecule.lineage_id, sep='\t', file=molecules_file)

    cells = compute_cells(corrected_molecules, args.min_length)
    logger.info(f'Detected {len(cells)} cells')
    write_cells(output_dir / 'cells.txt', cells)

    cells = filter_cells(cells, corrected_molecules, args.keep_single_reads)
    write_cells(output_dir / 'cells_filtered.txt', cells)

    lineages = compute_lineages(cells)
    logger.info(f'Detected {len(lineages)} lineages')
    lineage_sizes = Counter(len(cells) for cells in lineages.values())
    logger.info('Lineage size histogram (size: count): %s',
        ', '.join(f'{k}: {v}' for k, v in lineage_sizes.items()))
    with open(output_dir / 'lineages.txt', 'w') as f:
        print(
            '#Each output line corresponds to one barcode group (clone) and has '
            'the following style: Barcode\t:\tCellID1\tCellID2...\n'
            '# dash (-) = barcode base outside of read, '
            '0 = deletion in barcode sequence (position unknown)', file=f)

        for lineage_id in sorted(lineages):
            cells = lineages[lineage_id]
            row = [lineage_id, ':']
            for cell in cells:
                row.append(cell.cell_id)
            print(*row, sep='\t', file=f)

    # 2. Groups cells with same barcodes that most likely stem from one clone. Outputs a file
    #    with all clones and cellIDs belonging to each clone

    groups_dict = defaultdict(list)
    for cell in cells:
        # forms groups of cells with same lineage_id
        for lineage_id, count in cell.lineage_id_counts.items():
            groups_dict[lineage_id].append((cell.cell_id, count))

    logger.info(f'Detected {len(groups_dict)} unique cell groups')
    # in groups.txt all barcodes and their corresponding cellIDs can be found
    with open(output_dir / 'groups.txt', 'w') as groups_file:
        print(
            '#Each output line corresponds to one barcode group (clone) and has '
            'the following style: Barcode\t:\tCellID1\tbarcode-count1\tCellID2\tbarcode-count2...\n'
            '# dash (-) = barcode base outside of read, '
            '0 = deletion in barcode sequence (position unknown)', file=groups_file)

        for lineage_id in sorted(groups_dict):
            row = [lineage_id, ':']
            for cell_id, barcode_count in groups_dict[lineage_id]:
                row.extend([cell_id, barcode_count])
            print(*row, sep='\t', file=groups_file)

    # Create a loom file if requested
    if args.loom:
        write_loom(cells, input_dir, output_dir, lineage_id_length=args.end - args.start + 1)

    logger.info('Run completed!')


if __name__ == '__main__':
    main()

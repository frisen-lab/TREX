"""
Extract and filter random barcodes from single-cell sequencing data
"""
import sys
import math
import argparse
import operator
import shutil
import warnings
import logging
import subprocess
from io import StringIO
from pathlib import Path
from collections import Counter, defaultdict, OrderedDict
from typing import Set, List, Dict, NamedTuple, Iterable, Callable
from pkg_resources import get_distribution, DistributionNotFound

from alignlib import hamming_distance
import numpy as np
import pandas as pd
import pysam
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', 'Conversion of the second argument of issubdtype')
    import loompy

from .cellranger import make_cellranger_outs, CellRangerError
from .utils import NiceFormatter


try:
    __version__ = get_distribution('braintrace').version
except DistributionNotFound:
    # package is not installed
    pass

__author__ = 'leonie.von.berlin@ki.se'

logger = logging.getLogger(__name__)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    parser.add_argument('--genome-name', metavar='NAME',
        help='Name of the genome as indicated in Cell Ranger count run with the flag --genome. '
             'Default: Auto-detected',
        default=None)
    parser.add_argument('--chromosome', '--chr',
        help='Barcode chromosome name. Default: Last chromosome in the BAM file',
        default=None)
    parser.add_argument('--output', '-o', '--name', '-n', metavar='DIRECTORY', type=Path,
        help='name of the run and directory created by program. Default: %(default)s',
        default=Path('lineage_run'))
    parser.add_argument('--delete', action='store_true', help='Delete output directory if it exists')
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
    parser.add_argument('--amplicon', '-a', metavar='DIRECTORY',
        help='Path to Cell Ranger "outs" directory containing sequencing of the EGFP-barcode amplicon library. When combining' 
        'Cell Ranger runs indicate paths separated by comma and in same order as paths for transcriptome data. Example "path1,path2,path3"',
        default=None)
    parser.add_argument('--filter-cellids', '-f', metavar='DIRECTORY', type=Path,
        help='Path to a .csv file containing cell IDs to keep in the analysis. This flag enables to remove cells e.g. doublets',
        default=None)
    parser.add_argument('--keep-single-reads', action='store_true', default=False,
        help='Keep barcodes supported by only a single read. Default: Discard them')
    parser.add_argument('-l', '--loom',
        help='If given, create loom-file from Cell Ranger and barcode data. '
            'File will have the same name as the run',
        action='store_true')
    parser.add_argument('--restrict', metavar='FILE',
        help='Restrict analysis to the cell IDs listed in FILE')
    parser.add_argument('--highlight',
        help='Highlight cell IDs listed in FILE in the lineage graph')
    parser.add_argument('--cellid-suffix',
        help='Add suffixes to cell IDs to merge different Cell Ranger runs. Suffixes should be separated by comma and have the same order'
        'as the given Cell Ranger directory paths. Example: "_1,_2,_3"',
        default=None)
    parser.add_argument('--umi-matrix', default=False, action='store_true',
        help='Creates a umi count matrix with cells as columns and clone IDs as rows')
    parser.add_argument('--plot', dest='plot', default=False, action='store_true',
        help='Plot the lineage graph')
    parser.add_argument('path', metavar='DIRECTORY', type=Path,
        help='Path to Cell Ranger "outs" directory. To combine several runs, please separate paths by comma. Example: "path1,path2,path3".'
        'Do not forget to indicate cell IDs suffixes to separate cell IDs from differen runs with the --cellid-suffix flag')
    return parser.parse_args()


class Graph:
    """Undirected graph that can find connected components"""

    def __init__(self, nodes):
        # values are lists of adjacent nodes
        self._nodes = OrderedDict()
        for node in nodes:
            self._nodes[node] = []

    def add_edge(self, node1, node2):
        self._nodes[node1].append(node2)
        self._nodes[node2].append(node1)

    def remove_edge(self, node1, node2):
        self._nodes[node1].remove(node2)
        self._nodes[node2].remove(node1)

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

            subgraph = Graph(component)
            for n in component:
                subgraph._nodes[n] = self._nodes[n].copy()
            components.append(subgraph)
        return components

    def nodes(self):
        """Return all nodes as a list"""
        return list(self._nodes)

    def edges(self):
        """Yield all edges as pair (node1, node2)"""
        seen = set()
        for node1 in self._nodes:
            seen.add(node1)
            for node2 in self._nodes[node1]:
                if node2 in seen:
                    continue
                yield node1, node2

    def neighbors(self, node):
        """Return a list of all neighbors of a node"""
        return self._nodes[node]


def kmers(s: str, k: int):
    """
    Yield all overlapping k-kmers of s.

    >>> list(kmers('hello', 3))
    ['hel', 'ell', 'llo']
    """
    for i in range(len(s) - k + 1):
        yield s[i:i+k]


def cluster_sequences(
        sequences: List[str], is_similar: Callable[[str, str], bool], k: int = 6) -> List[List[str]]:
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
    return [g.nodes() for g in graph.connected_components()]


class Read(NamedTuple):
    umi: str
    cell_id: str
    clone_id: str


class Molecule(NamedTuple):
    umi: str
    cell_id: str
    clone_id: str
    read_count: int


class Cell(NamedTuple):
    cell_id: str
    clone_id_counts: Dict[str, int]

    def __hash__(self):
        return hash(self.cell_id)


class CellSet:
    def __init__(self, cells: List[Cell]):
        self.cells = cells
        self.cell_ids = tuple(sorted(c.cell_id for c in cells))
        self.clone_id_counts = sum((Counter(c.clone_id_counts) for c in cells), Counter())
        self.n = len(cells)
        self.cell_id = 'M-' + min(self.cell_ids)
        self._hash = hash(self.cell_ids)

    def __repr__(self):
        return f'CellSet(cells={self.cells!r})'

    def __hash__(self):
        return self._hash


def detect_clone_id_location(alignment_file, reference_name):
    """
    Detect where the clone id is located on the reference by inspecting the alignments.

    Return (clone_id_start, clone_id_end)
    """
    # Look for reference positions at which reads are soft-clipped at their 3' end
    starts = Counter()
    reference_length = alignment_file.get_reference_length(reference_name)
    for alignment in alignment_file.fetch(reference_name):
        clip_right = alignment.query_length - alignment.query_alignment_end
        if clip_right >= 5:
            starts[alignment.reference_end] += 1

    for clone_id_start, freq in starts.most_common(5):
        # Soft-clipping at the 5' end cannot be used to find the clone id end when
        # the clone id region is too far at the 3' end of the contig. Instead,
        # look at pileups and check base frequencies (over the clone id, bases should
        # be roughly uniformly distributed).
        if clone_id_start >= reference_length:
            # The most common reference position that is soft clipped is often the 3' end
            # of the contig. Skip that.
            continue
        clone_id_end = clone_id_start
        for column in alignment_file.pileup(reference_name, start=clone_id_start):
            if column.reference_pos < clone_id_start:
                # See pileup() documentation
                continue
            bases = [p.alignment.query_sequence[p.query_position] for p in column.pileups if p.query_position is not None]
            counter = Counter(bases)
            # Check whether one base dominates
            if counter.most_common()[0][1] / len(bases) > 0.95:
                # We appear to have found the end of the clone id
                clone_id_end = column.reference_pos
                break
        if clone_id_end - clone_id_start >= 5:
            # Good enough
            return (clone_id_start, clone_id_end)
    raise ValueError(f'Could not detect clone id location on chromosome {reference_name}')


def read_bam(bam_path: Path, output_dir: Path, allowed_cell_ids, chr_name, clone_id_start=None, clone_id_end=None, file_name_suffix="_entries", cellid_suffix=None):
    """
    bam_path -- path to input BAM file
    output_dir -- path to an output directory into which a BAM file is written that contais all
        reads on the chromosome that have the required tags.
    """
    cache_hit = cache_miss = 0
    with pysam.AlignmentFile(bam_path) as alignment_file:
        if chr_name is None:
            chr_name = alignment_file.references[-1]

        # TODO move out of here
        if clone_id_start is None or clone_id_end is None:
            if clone_id_start is not None or clone_id_end is not None:
                raise ValueError('Either both or none of clone id start and end must be provided')
            clone_id_start, clone_id_end = detect_clone_id_location(alignment_file, chr_name)
        logger.info(f"Reading clone ids from {chr_name}:{clone_id_start + 1}-{clone_id_end} "
            f"in {bam_path}")
        if clone_id_end - clone_id_start < 10:
            raise ValueError('Auto-detected clone id too short, something is wrong')
        output_bam_path = output_dir / (chr_name + file_name_suffix + '.bam')

        with pysam.AlignmentFile(output_bam_path, 'wb', template=alignment_file) as out_bam:
            # Fetches those reads aligning to the artifical, clone-id-containing chromosome
            reads = []
            unknown_ids = no_cell_id = no_umi = 0
            prev_start = None
            for read in alignment_file.fetch(chr_name, max(0, clone_id_start - 10), clone_id_end + 10):
                # Skip reads without cellID or UMI
                if not read.has_tag('CB') or not read.has_tag('UB'):
                    if not read.has_tag('CB'):
                        no_cell_id += 1
                    if not read.has_tag('UB'):
                        no_umi += 1
                    continue
                # Filters out reads that have not approved cellIDs
                cell_id = read.get_tag('CB')
                if cell_id not in allowed_cell_ids:
                    unknown_ids += 1
                    continue
                # Gives cellIDs an experiment-specific suffix
                if cellid_suffix is not None:
                    if cell_id.endswith("-1"):
                        cell_id = cell_id[:-2]
                    cell_id += cellid_suffix

                # Use a cache of clone id extraction results. This helps when there are
                # many identical reads (amplicon data)
                if prev_start != read.reference_start:
                    clone_id_cache = dict()
                    prev_start = read.reference_start
                cache_key = (read.reference_start, read.query_sequence)
                try:
                    clone_id = clone_id_cache[cache_key]
                    cache_hit += 1
                except KeyError:
                    cache_miss += 1
                    clone_id = extract_clone_id(clone_id_start, clone_id_end, read)
                    clone_id_cache[cache_key] = clone_id
                if clone_id is None:
                    # Read does not cover the clone id
                    continue
                reads.append(Read(cell_id=cell_id, umi=read.get_tag('UB'), clone_id=clone_id))

                # Write the passing alignments to a separate file
                out_bam.write(read)

            logger.info(f'Skipped {unknown_ids} reads with unrecognized cell ids '
                        f'(and {no_umi+no_cell_id} without UMI or cell id)')
    logger.info(f"Cache hits: {cache_hit}. Cache misses: {cache_miss}")
    sorted_reads = sorted(reads, key=lambda read: (read.umi, read.cell_id, read.clone_id))
    assert len(sorted_reads) == 0 or len(sorted_reads[0].clone_id) == clone_id_end - clone_id_start
    return sorted_reads


def extract_clone_id(clone_id_start, clone_id_end, read):
    query_align_end = read.query_alignment_end
    query_align_start = read.query_alignment_start
    query_sequence = read.query_sequence
    # Extract clone id
    clone_id = ['-'] * (clone_id_end - clone_id_start)
    bases = 0
    for query_pos, ref_pos in read.get_aligned_pairs():
        # Replace soft-clipping with an ungapped alignment extending into the
        # soft-clipped region, assuming the clipping occurred because the clone id
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

        if ref_pos is not None and clone_id_start <= ref_pos < clone_id_end:
            if query_pos is None:
                # Deletion or intron skip
                query_base = '0'
            else:
                # Match or mismatch
                query_base = query_sequence[query_pos]
                bases += 1
            clone_id[ref_pos - clone_id_start] = query_base

    if bases == 0:
        return None
    else:
        return "".join(clone_id)


def read_allowed_cellids(cellids_filtered):
    """
    Read a user-provided list of allowed cellIDs from Seurat like this:

    AAACCTGAGCGACGTA

    OR:

    AAACCTGCATACTCTT_1
    """
    allowed_ids = []
    filtered_df = pd.read_csv(Path(cellids_filtered), sep=",", index_col=0)
    for line in filtered_df.iloc[:, 0]:
        allowed_ids.append(line.split("_")[0] + "-1")
    logger.info(f'Restricting analysis to {len(allowed_ids)} allowed cells')
    return set(allowed_ids)


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

    - forms consensus sequence of all clone ids of one group,
    """
    groups = defaultdict(list)
    for read in sorted_reads:
        groups[(read.umi, read.cell_id)].append(read.clone_id)

    molecules = []
    for (umi, cell_id), clone_ids in groups.items():
        clone_id_consensus = compute_consensus(clone_ids)
        molecules.append(
            Molecule(cell_id=cell_id, umi=umi, clone_id=clone_id_consensus,
                read_count=len(clone_ids)))

    sorted_molecules = sorted(molecules, key=lambda mol: (mol.cell_id, mol.clone_id, mol.umi))

    return sorted_molecules


def correct_clone_ids(
        molecules: List[Molecule], max_hamming: int, min_overlap: int = 20) -> List[Molecule]:
    """
    Attempt to correct sequencing errors in the clone id sequences of all molecules
    """
    # Obtain all clone ids (including those with '-' and '0')
    clone_ids = [m.clone_id for m in molecules]

    # Count the full-length clone ids
    clone_id_counts = Counter(clone_ids)

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

    clusters = cluster_sequences(list(set(clone_ids)), is_similar=is_similar, k=7)

    # Map non-singleton clone ids to a cluster representative
    clone_id_map = dict()
    for cluster in clusters:
        if len(cluster) > 1:
            # Pick most frequent clone id as representative
            representative = max(cluster, key=lambda bc: (clone_id_counts[bc], bc))
            for clone_id in cluster:
                clone_id_map[clone_id] = representative

    # Create a new list of molecules in which the clone ids have been replaced
    # by their representatives
    new_molecules = []
    for molecule in molecules:
        clone_id = clone_id_map.get(molecule.clone_id, molecule.clone_id)
        new_molecules.append(molecule._replace(clone_id=clone_id))
    return new_molecules


def compute_cells(sorted_molecules: List[Molecule], minimum_clone_id_length: int) -> List[Cell]:
    """
    Group molecules into cells.
    """
    # 1. Forms groups of molecules (with set clone id minimum length) that have identical cellIDs
    #    => belong to one cell,
    # 2. counts number of appearances of each clone id in each group,

    cell_id_groups = defaultdict(list)
    for molecule in sorted_molecules:
        clone_id = molecule.clone_id
        pure_li = clone_id.strip('-')
        # TODO may not work as intended (strip only removes prefixes and suffixes)
        pure_bc0 = clone_id.strip('0')
        if len(pure_li) >= minimum_clone_id_length and len(pure_bc0) >= minimum_clone_id_length:
            cell_id_groups[molecule.cell_id].append(molecule)

    cells = []
    for cell_id, molecules in cell_id_groups.items():
        clone_ids = [molecule.clone_id for molecule in molecules]
        clone_id_counts = OrderedDict(sorted(Counter(clone_ids).most_common(),
            key=lambda x: x[0].count('-')))
        cells.append(Cell(cell_id=cell_id, clone_id_counts=clone_id_counts))
    return cells


def filter_cells(
        cells: Iterable[Cell], molecules: Iterable[Molecule],
        keep_single_reads: bool = False) -> List[Cell]:
    """
    Filter clone ids according to two criteria:

    - Clone ids that have only a count of one and can be found in another cell are most
      likely results of contamination and are removed,
    - If keep_single_reads is False, clone ids that have only a count of one and are also only based
      on one read are also removed
    """
    overall_clone_id_counts: Dict[str, int] = Counter()
    for cell in cells:
        overall_clone_id_counts.update(cell.clone_id_counts)

    single_read_clone_ids = set()
    for molecule in molecules:
        if molecule.read_count == 1:
            single_read_clone_ids.add(molecule.clone_id)
    logger.info(f"Found {len(single_read_clone_ids)} single-read clone ids")

    # filters out clone ids with a count of one that appear in another cell
    new_cells = []
    for cell in cells:
        clone_id_counts = cell.clone_id_counts.copy()
        for clone_id, count in cell.clone_id_counts.items():
            if count > 1:
                # This clone id occurs more than once in this cell - keep it
                continue
            if overall_clone_id_counts[clone_id] > 1:
                # This clone id occurs also in other cells - remove it
                del clone_id_counts[clone_id]
            elif clone_id in single_read_clone_ids and not keep_single_reads:
                del clone_id_counts[clone_id]
        if clone_id_counts:
            new_cells.append(Cell(cell_id=cell.cell_id, clone_id_counts=clone_id_counts))
    return new_cells


class LineageGraph:
    def __init__(self, cells: List[Cell]):
        self._cells = cells
        self._graph = self._make_cell_graph()

    def _make_cell_graph(self):
        """
        Create graph of cells; add edges between cells that
        share at least one barcode. Return created graph.
        """
        cells = [cell for cell in self._cells if cell.clone_id_counts]
        graph = Graph(cells)
        for i in range(len(cells)):
            for j in range(i + 1, len(cells)):
                if set(cells[i].clone_id_counts) & set(cells[j].clone_id_counts):
                    # Cell i and j share a clone id
                    graph.add_edge(cells[i], cells[j])
        return graph

    def _make_barcode_graph(self):
        clone_id_counts: Dict[str, int] = Counter()
        for cell in self._cells:
            clone_id_counts.update(cell.clone_id_counts)
        all_clone_ids = list(clone_id_counts)

        # Create a graph of barcodes; add an edge for barcodes occuring in the same cell
        graph = Graph(all_clone_ids)
        for cell in self._cells:
            barcodes = list(cell.clone_id_counts)
            for other in barcodes[1:]:
                graph.add_edge(barcodes[0], other)
        return graph

    def lineages(self) -> Dict[str, List[Cell]]:
        """
        Compute lineages. Return a dict that maps a representative clone id to a list of cells.
        """
        clusters = [g.nodes() for g in self._graph.connected_components()]

        def most_abundant_clone_id(cells: List[Cell]):
            counts = Counter()
            for cell in cells:
                counts.update(cell.clone_id_counts)
            return max(counts, key=lambda k: (counts[k], k))

        return {most_abundant_clone_id(cells): cells for cells in clusters}

    def dot(self, highlight=None):
        s = StringIO()
        print('graph g {', file=s)
        for node1, node2 in self._graph.edges():
            print(f'"{node1.cell_id}" -- "{node2.cell_id}"', file=s)
        print('}', file=s)
        return s.getvalue()

    @property
    def graph(self):
        return self._graph


class CompressedLineageGraph:
    def __init__(self, cells: List[Cell]):
        self._cells = self._compress_cells(cells)
        self._graph = self._make_cell_graph()

    @staticmethod
    def _compress_cells(cells):
        cell_lists = defaultdict(list)
        for cell in cells:
            clone_ids = tuple(cell.clone_id_counts)
            cell_lists[clone_ids].append(cell)

        cell_sets = []
        for cells in cell_lists.values():
            cell_sets.append(CellSet(cells))
        return cell_sets

    def _make_cell_graph(self):
        """
        Create graph of cells; add edges between cells that
        share at least one clone id. Return created graph.
        """
        cells = [cell for cell in self._cells if cell.clone_id_counts]
        graph = Graph(cells)
        for i in range(len(cells)):
            for j in range(i + 1, len(cells)):
                if set(cells[i].clone_id_counts) & set(cells[j].clone_id_counts):
                    # Cell i and j share a clone id
                    graph.add_edge(cells[i], cells[j])
        return graph

    def bridges(self):
        """Find edges that appear to incorrectly bridge two unrelated subclusters."""
        # A bridge as defined here is simply an edge between two nodes that
        # have a non-empty set of common neighbors. We do not check for actual
        # connectivity at the moment.
        bridges = []
        for node1, node2 in self._graph.edges():
            neighbors1 = self._graph.neighbors(node1)
            neighbors2 = self._graph.neighbors(node2)
            common_neighbors = set(neighbors1) & set(neighbors2)
            if not common_neighbors and (len(neighbors1) > 1 or len(neighbors2) > 1):
                bridges.append((node1, node2))
        return bridges

    def remove_edges(self, edges):
        for node1, node2 in edges:
            self._graph.remove_edge(node1, node2)

    @staticmethod
    def _expand_cell_sets(cell_sets: List[CellSet]) -> List[Cell]:
        """Expand a list of CellSets into a list of Cells"""
        cells = []
        for cell_set in cell_sets:
            cells.extend(cell_set.cells)
        return cells

    def write_lineages(self, path):
        lineages = self.lineages()
        with open(path, 'w') as f:
            print(
                '#Each output line corresponds to one barcode group (clone) and has '
                'the following style: Barcode\t:\tCellID1\tCellID2...\n'
                '# dash (-) = barcode base outside of read, '
                '0 = deletion in barcode sequence (position unknown)', file=f)

            for clone_id in sorted(lineages):
                cells = sorted(lineages[clone_id])
                row = [clone_id, ':']
                for cell in cells:
                    row.append(cell.cell_id)
                print(*row, sep='\t', file=f)
        return lineages

    def lineages(self) -> Dict[str, List[Cell]]:
        """
        Compute lineages. Return a dict that maps a representative clone id to a list of cells.
        """
        compressed_clusters = [g.nodes() for g in self._graph.connected_components()]
        # Expand the CellSet instances into cells
        clusters = [self._expand_cell_sets(cluster) for cluster in compressed_clusters]

        def most_abundant_clone_id(cells: List[Cell]):
            counts = Counter()
            for cell in cells:
                counts.update(cell.clone_id_counts)
            return max(counts, key=lambda k: (counts[k], k))

        return {most_abundant_clone_id(cells): cells for cells in clusters}

    def plot(self, path, highlight=None):
        graphviz_path = path.with_suffix(".gv")
        with open(graphviz_path, "w") as f:
            print(self.dot(highlight), file=f)
        pdf_path = str(path.with_suffix(".pdf"))
        subprocess.run(["sfdp", "-Tpdf", "-o", pdf_path, graphviz_path], check=True)

    def dot(self, highlight=None):
        if highlight is not None:
            highlight = set(highlight)
        max_width = 10
        edge_scaling = (max_width - 1) / math.log(
            max((node1.n * node2.n for node1, node2 in self._graph.edges()), default=math.exp(1)))
        node_scaling = (max_width - 1) / math.log(max(node.n for node in self._graph.nodes()))
        s = StringIO()
        print('graph g {', file=s)
        # Using overlap=false would be nice here, but that does not work with some Graphviz builds
        print('  graph [outputorder=edgesfirst];', file=s)
        print('  edge [color=blue];', file=s)
        print('  node [style=filled, fillcolor=white, fontname="Roboto"];', file=s)
        for node in self._graph.nodes():
            if self._graph.neighbors(node):
                width = int(1 + node_scaling * math.log(node.n))
                intersection = set(node.cell_ids) & highlight
                hl = ',fillcolor=yellow' if intersection else ''
                hl_label = f' ({len(intersection)})' if intersection else ''
                print(
                    f'  "{node.cell_id}" [penwidth={width}{hl},label="{node.cell_id}\\n{node.n}{hl_label}"];',
                    file=s)
        for node1, node2 in self._graph.edges():
            width = int(1 + edge_scaling * math.log(node1.n * node2.n))
            neighbors1 = self._graph.neighbors(node1)
            neighbors2 = self._graph.neighbors(node2)
            common_neighbors = set(neighbors1) & set(neighbors2)
            bridge = ''
            if (len(neighbors1) > 1 or len(neighbors2) > 1) and not common_neighbors:
                bridge = ', style=dashed, color=red'
                width = 2
            print(f'  "{node1.cell_id}" -- "{node2.cell_id}" [penwidth={width}{bridge}];', file=s)

        print('}', file=s)
        return s.getvalue()

    def components_txt(self, highlight=None):
        s = StringIO()
        print('# Lineage graph components (only incomplete/density<1)', file=s)
        n_complete = 0
        for subgraph in self.graph.connected_components():
            cells = sorted(self._expand_cell_sets(subgraph.nodes()), key=lambda c: c.cell_id)
            n_nodes = len(cells)
            n_edges = sum(n1.n * n2.n for n1, n2 in subgraph.edges())
            n_edges += sum(node.n * (node.n - 1) // 2 for node in subgraph.nodes())
            possible_edges = n_nodes * (n_nodes - 1) // 2
            if n_edges == possible_edges:
                n_complete += 1
                continue
            density = n_edges / possible_edges
            print(f'## {n_nodes} nodes, {n_edges} edges, density {density:.3f}', file=s)
            counter = Counter()
            for cell in cells:
                if highlight is not None and cell.cell_id in highlight:
                    highlighting = '+'
                else:
                    highlighting = ''
                print(cell.cell_id, highlighting, *sorted(cell.clone_id_counts.keys()), sep='\t', file=s)
                counter.update(cell.clone_id_counts.keys())
        print(f'# {n_complete} complete components', file=s)
        return s.getvalue()

    @property
    def graph(self):
        return self._graph


def write_loom(cells: List[Cell], cellranger_outs, output_dir, clone_id_length, top_n=6):
    """
    Create a loom file from a Cell Ranger sample directory and augment it with information about
    the most abundant clone id and their counts.
    """
    # For each cell, collect the most abundant clone ids and their counts
    # Maps cell_id to a list of (clone_id, count) pairs that represent the most abundant clone ids.
    most_abundant = dict()
    for cell in cells:
        if not cell.clone_id_counts:
            continue
        counts = sorted(cell.clone_id_counts.items(), key=operator.itemgetter(1))
        counts.reverse()
        counts = counts[:top_n]
        most_abundant[cell.cell_id] = counts

    loompy.create_from_cellranger(cellranger_outs.sample_dir, outdir=output_dir)
    # create_from_cellranger() does not tell us the name of the created file,
    # so we need to re-derive it from the sample name.
    sample_name = cellranger_outs.sample_dir.name
    loom_path = output_dir / (sample_name + '.loom')

    with loompy.connect(loom_path) as ds:
        # Cell ids in the loom file are prefixed by the sample name and a ':'. Remove that prefix.
        loom_cell_ids = [cell_id[len(sample_name)+1:] for cell_id in ds.ca.CellID]

        # Transform clone ids and count data
        # brings clone id data into correct format for loom file.
        # Array must have same shape as all_cellIDs
        clone_id_lists = [[] for _ in range(top_n)]
        count_lists = [[] for _ in range(top_n)]
        for cell_id in loom_cell_ids:
            clone_id_counts = most_abundant.get(cell_id, [])
            # Fill up to a constant length
            while len(clone_id_counts) < top_n:
                clone_id_counts.append(('-', 0))

            for i, (clone_id, count) in enumerate(clone_id_counts):
                clone_id_lists[i].append(clone_id)
                count_lists[i].append(count)

        # Add clone id and count information to loom file
        for i in range(top_n):
            ds.ca[f'linBarcode_{i+1}'] = np.array(clone_id_lists[i], dtype='S%r' % clone_id_length)
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
            sorted_clone_ids = sorted(cell.clone_id_counts, key=lambda x: cell.clone_id_counts[x], reverse=True)
            if not sorted_clone_ids:
                continue
            for clone_id in sorted_clone_ids:
                row.extend([clone_id, cell.clone_id_counts[clone_id]])
            print(*row, sep='\t', file=f)


def write_umimatrix(output_dir: Path, cells: List[Cell]):
    """Create a UMI-count matrix with cells as columns and clone ids as rows"""
    clone_ids = set()
    for cell in cells:
        clone_ids.update(clone_id for clone_id in cell.clone_id_counts)
    clone_ids = sorted(clone_ids)
    all_clone_id_counts = [cell.clone_id_counts for cell in cells]
    with open(output_dir / "umi_count_matrix.csv", "w") as f:
        f.write(",")
        f.write(",".join(cell.cell_id for cell in cells))
        f.write("\n")
        for clone_id in clone_ids:
            f.write(clone_id)
            f.write(",")
            values = [lic.get(clone_id, 0) for lic in all_clone_id_counts]
            f.write(",".join(str(v) for v in values))
            f.write("\n")


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
    output_dir = args.output
    setup_logging(debug=False)

    # PART I + II: Clone id extraction and reads construction

    # 1. Extracts reads aligning to clone id chromosome,
    # 2. extracts clone ids, UMIs and cellIDs from reads,
    # 3. outputs UMI-sorted reads with barcodes

    try:
        make_output_dir(output_dir, args.delete)
    except FileExistsError:
        logger.error(f'Output directory "{output_dir}" already exists '
                     '(use --delete to force deleting an existing output directory)')
        sys.exit(1)

    add_file_logging(output_dir / 'log.txt')
    logger.info('Command line arguments: %s', ' '.join(sys.argv[1:]))

    restrict_cell_ids = None
    if args.restrict:
        with open(args.restrict) as f:
            restrict_cell_ids = [line.strip() for line in f]
    allowed_cell_ids = None
    if args.filter_cellids:
        allowed_cell_ids = read_allowed_cellids(args.filter_cellids)
    transcriptome_inputs = str(args.path).split(",")
    if args.cellid_suffix:
        cellid_suffixes = args.cellid_suffix.split(",")
    else:
        if len(transcriptome_inputs) == 1:
            cellid_suffixes = [None]  # Do not modify suffixes
        else:
            cellid_suffixes = ["_{}".format(i) for i in range(1, len(transcriptome_inputs) + 1)]
            logger.info("Using these cellid suffixes: %s", ", ".join(cellid_suffixes))
    if len(cellid_suffixes) != len(transcriptome_inputs):
        logger.error("The number of cellid suffixs (--cellid-suffix) must match the number of "
            "provided transcriptome datasets")
        sys.exit(1)
    if args.amplicon:
        amplicon_inputs = [Path(a) for a in args.amplicon.split(",")]
        if len(transcriptome_inputs) != len(amplicon_inputs):
            logger.error("As many amplicon as transcriptome datasets must be provided")
            sys.exit(1)
    else:
        amplicon_inputs = []

    def read_one_dataset(path, suffix, file_name_suffix):
        outs_dir = make_cellranger_outs(path, args.genome_name)
        if allowed_cell_ids:
            cell_ids = allowed_cell_ids
        else:
            cell_ids = outs_dir.cellids()
        logger.info(f'Found {len(cell_ids)} cell ids in the barcodes.tsv file')

        return read_bam(
            outs_dir.bam, output_dir, cell_ids, args.chromosome,
            args.start - 1 if args.start is not None else None, args.end,
            file_name_suffix=file_name_suffix, cellid_suffix=suffix)

    # Extracts reads from  and amplicon clone id chromosome from amplicon sequencing data
    # Combines reads from amplicon dataset with reads from transcriptome dataset for
    # clone id, cellID and UMI extraction
    sorted_reads = list()
    try:
        for path, suffix in zip(transcriptome_inputs, cellid_suffixes):
            sorted_reads.extend(read_one_dataset(path, suffix, "_entries"))
        for path, suffix in zip(amplicon_inputs, cellid_suffixes):
            sorted_reads.extend(read_one_dataset(path, suffix, "_amp_entries"))
    except CellRangerError as e:
        logger.error("%s", e)
        sys.exit(1)
    sorted_reads.sort(key=lambda read: (read.umi, read.cell_id, read.clone_id))

    clone_ids = [
        r.clone_id for r in sorted_reads if '-' not in r.clone_id and '0' not in r.clone_id]
    logger.info(f'Read {len(sorted_reads)} reads containing (parts of) the barcode '
        f'({len(clone_ids)} full barcodes, {len(set(clone_ids))} unique)')

    write_reads(output_dir / "reads.txt", sorted_reads)

    # Part III: Molecule construction

    molecules = compute_molecules(sorted_reads)
    clone_ids = [
        m.clone_id for m in molecules if '-' not in m.clone_id and '0' not in m.clone_id]
    logger.info(f'Detected {len(molecules)} molecules ({len(clone_ids)} full clone ids, '
        f'{len(set(clone_ids))} unique)')

    write_molecules(output_dir / 'molecules.txt', molecules)

    corrected_molecules = correct_clone_ids(molecules, args.max_hamming, args.min_length)
    clone_ids = [m.clone_id for m in corrected_molecules
        if '-' not in m.clone_id and '0' not in m.clone_id]
    logger.info(f'After clone id correction, {len(set(clone_ids))} unique clone ids remain')

    write_molecules(output_dir / 'molecules_corrected.txt', corrected_molecules)

    cells = compute_cells(corrected_molecules, args.min_length)
    logger.info(f'Detected {len(cells)} cells')
    write_cells(output_dir / 'cells.txt', cells)

    cells = filter_cells(cells, corrected_molecules, args.keep_single_reads)
    logger.info(f'{len(cells)} filtered cells remain')
    write_cells(output_dir / 'cells_filtered.txt', cells)

    if args.umi_matrix:
        logger.info(f"Writing UMI matrix")
        write_umimatrix(output_dir, cells)

    if restrict_cell_ids is not None:
        restrict_cell_ids = set(restrict_cell_ids)
        cells = [cell for cell in cells if cell.cell_id in restrict_cell_ids]
        logger.info(f'Restricting to {len(cells)} cells')

    lineage_graph = CompressedLineageGraph(cells)
    highlight_cell_ids = []
    if args.highlight:
        with open(args.highlight) as f:
            highlight_cell_ids = [line.strip() for line in f]

    with open(output_dir / 'components.txt', 'w') as components_file:
        print(lineage_graph.components_txt(highlight_cell_ids), file=components_file, end='')
    if args.plot:
        logger.info('Plotting compressed lineage graph')
        lineage_graph.plot(output_dir / 'graph', highlight_cell_ids)

    bridges = lineage_graph.bridges()
    logger.info(f'Removing {len(bridges)} bridges from the graph')
    lineage_graph.remove_edges(bridges)
    with open(output_dir / 'components_corrected.txt', 'w') as components_file:
        print(lineage_graph.components_txt(highlight_cell_ids), file=components_file, end='')

    if args.plot:
        logger.info('Plotting corrected lineage graph')
        lineage_graph.plot(output_dir / 'graph_corrected', highlight_cell_ids)

    lineages = lineage_graph.write_lineages(output_dir / 'lineages.txt')
    logger.info(f'Detected {len(lineages)} lineages')
    lineage_sizes = Counter(len(cells) for cells in lineages.values())
    logger.info('Lineage size histogram (size: count): %s',
        ', '.join(f'{k}: {v}' for k, v in lineage_sizes.items()))

    if args.loom:
        if len(transcriptome_inputs) > 1:
            logger.warning("Writing a loom file only for the first transcriptome dataset")
        outs_dir = make_cellranger_outs(transcriptome_inputs[0])
        write_loom(cells, outs_dir, output_dir, clone_id_length=args.end - args.start + 1)

    logger.info('Run completed!')


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


def write_reads(path, reads):
    with open(path, 'w') as f:
        print(
            '#Each output line corresponds to one read and has the following style: '
            'CellID\tUMI\tBarcode\n'
            '# dash (-) = barcode base outside of read, '
            '0 = deletion in barcode sequence (position unknown)', file=f)
        for read in reads:
            print(read.cell_id, read.umi, read.clone_id, sep='\t', file=f)


def write_molecules(path, molecules):
    with open(path, 'w') as f:
        print(
            '#Each output line corresponds to one molecule and has the following style: '
            'CellID\tUMI\tBarcode\n'
            '# dash (-) = barcode base outside of read, '
            '0 = deletion in barcode sequence (position unknown)', file=f)
        for molecule in molecules:
            print(molecule.cell_id, molecule.umi, molecule.clone_id, sep='\t', file=f)


if __name__ == '__main__':
    main()

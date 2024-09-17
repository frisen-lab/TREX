"""
Clone computation. A clone is represented as a set of cells.
"""

from typing import List, Dict, Tuple
from collections import Counter, defaultdict
from io import StringIO
import math
import subprocess
import logging

from .graph import Graph
from .cell import Cell


logger = logging.getLogger(__name__)


class Clone:
    def __init__(self, cells: List[Cell]):
        self.cells = cells
        self.cell_ids = tuple(sorted(c.cell_id for c in cells))
        self.counts: Counter = sum((Counter(c.counts) for c in cells), Counter())
        self.n = len(cells)
        self.cell_id = "M-" + min(self.cell_ids)
        self._hash = hash(self.cell_ids)

    def __repr__(self):
        return f"Clone(cells={self.cells!r})"

    def __hash__(self):
        return self._hash


class CellGraph:
    """
    Graph of cells. An edge is drawn between two cells if the Jaccard similarity between
    their sets of cloneIDs is above the given threshold.
    """

    def __init__(self, cells: List[Cell], jaccard_threshold: float):
        self._jaccard_threshold = jaccard_threshold
        self._clones = self._precluster_cells(cells)
        self._graph = self._make_graph()

    @staticmethod
    def _precluster_cells(cells):
        """Put cells that have identical sets of cloneIDs into a clone"""
        cell_lists = defaultdict(list)
        for cell in cells:
            clone_ids = tuple(sorted(cell.counts))
            cell_lists[clone_ids].append(cell)

        clones = []
        for cells in cell_lists.values():
            clones.append(Clone(cells))
        return clones

    def _make_graph(self):
        """
        Create graph of clones; add edges between pre-clustered clones that
        share at least one cloneID. Return created graph.
        """
        clones = [clone for clone in self._clones if clone.counts]
        graph = Graph(clones)
        n = 0
        for i in range(len(clones)):
            for j in range(i + 1, len(clones)):
                if self._is_similar(clones[i], clones[j]):
                    n += 1
                    graph.add_edge(clones[i], clones[j])
        logger.debug(f"Added {n} edges to the clone graph")
        return graph

    @staticmethod
    def jaccard_index(a, b):
        if not a and not b:
            return 1
        return len(a & b) / len(a | b)

    def _is_similar(self, clone1, clone2):
        # TODO compute a weighted index using counts?
        a = set(clone1.counts)
        b = set(clone2.counts)
        index = self.jaccard_index(a, b)
        return index > self._jaccard_threshold

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

    def doublets(self):
        """Find cells that appear to incorrectly connect two unrelated subclusters"""
        cut_vertices = self._graph.local_cut_vertices()
        # Skip nodes that represent multiple cells
        return [node for node in cut_vertices if node.n == 1]

    def remove_edges(self, edges):
        for node1, node2 in edges:
            self._graph.remove_edge(node1, node2)

    def remove_nodes(self, nodes):
        for node in nodes:
            self._graph.remove_node(node)

    @staticmethod
    def _expand_clones(clones: List[Clone]) -> List[Cell]:
        """Expand a list of Clone instances into a list of Cells"""
        cells = []
        for clone in clones:
            cells.extend(clone.cells)
        return cells

    @staticmethod
    def write_clones(file, clones):
        print("clone_nr", "cell_id", sep="\t", file=file)
        for index, (clone_id, cells) in enumerate(sorted(clones), start=1):
            cells = sorted(cells)
            for cell in cells:
                print(index, cell.cell_id, sep="\t", file=file)

    @staticmethod
    def write_clone_sequences(file, clones):
        print("clone_nr", "clone_id", sep="\t", file=file)
        for index, (clone_id, cells) in enumerate(sorted(clones), start=1):
            print(index, clone_id, sep="\t", file=file)

    def clones(self) -> List[Tuple[str, List[Cell]]]:
        """
        Compute clones. Return a dict that maps a cloneID to a list of cells.
        """
        compressed_clusters = [g.nodes() for g in self._graph.connected_components()]
        # Expand the Clone instances into cells
        clusters = [self._expand_clones(cluster) for cluster in compressed_clusters]

        logger.debug(
            "CloneGraph.clones() called. %s connected components (clones) found",
            len(clusters),
        )

        def most_abundant_clone_id(cells: List[Cell]):
            counts: Counter = Counter()
            for cell in cells:
                counts.update(cell.counts)
            return max(counts, key=lambda k: (counts[k], k))

        return [(most_abundant_clone_id(cells), cells) for cells in clusters]

    def plot(self, path, highlight_cell_ids=None, highlight_doublets=None):
        graphviz_path = path.with_suffix(".gv")
        with open(graphviz_path, "w") as f:
            print(self.dot(highlight_cell_ids, highlight_doublets), file=f)
        pdf_path = str(path.with_suffix(".pdf"))
        subprocess.run(["sfdp", "-Tpdf", "-o", pdf_path, graphviz_path], check=True)

    def dot(self, highlight_cell_ids=None, highlight_doublets=None) -> str:
        highlight_cell_ids = (
            set(highlight_cell_ids) if highlight_cell_ids is not None else set()
        )
        highlight_doublets = (
            set(highlight_doublets) if highlight_doublets is not None else set()
        )
        max_width = 10
        try:
            edge_scaling = (max_width - 1) / math.log(
                max(
                    (node1.n * node2.n for node1, node2 in self._graph.edges()),
                    default=math.exp(1),
                )
            )
            node_scaling = (max_width - 1) / math.log(
                max(node.n for node in self._graph.nodes())
            )
        except ZeroDivisionError:
            edge_scaling = max_width - 1
            node_scaling = max_width - 1
        s = StringIO()
        print("graph g {", file=s)
        # Using overlap=false would be nice here,
        # but that does not work with some Graphviz builds
        print("  graph [outputorder=edgesfirst];", file=s)
        print("  edge [color=blue];", file=s)
        print('  node [style=filled, fillcolor=white, fontname="Roboto"];', file=s)
        for node in self._graph.nodes():
            if self._graph.neighbors(node):
                width = int(1 + node_scaling * math.log(node.n))
                intersection = set(node.cell_ids) & highlight_cell_ids
                if node in highlight_doublets:
                    hl = ",fillcolor=orange"
                    hl_label = " (doublet)"
                elif intersection:
                    hl = ",fillcolor=yellow"
                    hl_label = f" ({len(intersection)})"
                else:
                    hl = ""
                    hl_label = ""
                print(
                    f'  "{node.cell_id}" [penwidth={width}{hl},label="{node.cell_id}'
                    f'\\n{node.n}{hl_label}"];',
                    file=s,
                )
        for node1, node2 in self._graph.edges():
            width = int(1 + edge_scaling * math.log(node1.n * node2.n))
            neighbors1 = self._graph.neighbors(node1)
            neighbors2 = self._graph.neighbors(node2)
            common_neighbors = set(neighbors1) & set(neighbors2)
            bridge = ""
            if (len(neighbors1) > 1 or len(neighbors2) > 1) and not common_neighbors:
                bridge = ", style=dashed, color=red"
                width = 2
            print(
                f'  "{node1.cell_id}" -- "{node2.cell_id}" [penwidth={width}{bridge}];',
                file=s,
            )
        print("}", file=s)

        return s.getvalue()

    def components_txt(self, highlight=None):
        s = StringIO()
        print("# Clone graph components (only incomplete/density<1)", file=s)
        n_complete = 0
        for subgraph in self.graph.connected_components():
            cells = sorted(
                self._expand_clones(subgraph.nodes()), key=lambda c: c.cell_id
            )
            n_nodes = len(cells)
            n_edges = sum(n1.n * n2.n for n1, n2 in subgraph.edges())
            n_edges += sum(node.n * (node.n - 1) // 2 for node in subgraph.nodes())
            possible_edges = n_nodes * (n_nodes - 1) // 2
            if n_edges == possible_edges:
                n_complete += 1
                continue
            density = n_edges / possible_edges
            print(f"## {n_nodes} nodes, {n_edges} edges, density {density:.3f}", file=s)
            counter = Counter()
            for cell in cells:
                if highlight is not None and cell.cell_id in highlight:
                    highlighting = "+"
                else:
                    highlighting = ""
                clone_ids = sorted(cell.counts.keys())
                print(cell.cell_id, highlighting, *clone_ids, sep="\t", file=s)
                counter.update(cell.counts.keys())
        print(f"# {n_complete} complete components", file=s)
        return s.getvalue()

    @property
    def graph(self):
        return self._graph


# TODO this is unused
class UncompressedCellGraph:
    def __init__(self, cells: List[Cell]):
        self._cells = cells
        self._graph = self._make_cell_graph()

    def _make_cell_graph(self):
        """
        Create graph of cells; add edges between cells that
        share at least one barcode. Return created graph.
        """
        cells = [cell for cell in self._cells if cell.counts]
        graph = Graph(cells)
        for i in range(len(cells)):
            for j in range(i + 1, len(cells)):
                if set(cells[i].counts) & set(cells[j].counts):
                    # Cell i and j share a cloneID
                    graph.add_edge(cells[i], cells[j])
        return graph

    def _make_barcode_graph(self) -> Graph:
        counts: Dict[str, int] = Counter()
        for cell in self._cells:
            counts.update(cell.counts)
        all_clone_ids = list(counts)

        # Create a graph of barcodes; add an edge for barcodes occuring in the same cell
        graph = Graph(all_clone_ids)
        for cell in self._cells:
            barcodes = list(cell.counts)
            for other in barcodes[1:]:
                graph.add_edge(barcodes[0], other)
        return graph

    def clones(self) -> Dict[str, List[Cell]]:
        """
        Compute clones. Return a dict that maps a representative cloneID to a list of cells.
        """
        clusters = [g.nodes() for g in self._graph.connected_components()]

        def most_abundant_clone_id(cells: List[Cell]):
            counts: Counter = Counter()
            for cell in cells:
                counts.update(cell.counts)
            return max(counts, key=lambda k: (counts[k], k))

        return {most_abundant_clone_id(cells): cells for cells in clusters}

    def dot(self, highlight=None):
        s = StringIO()
        print("graph g {", file=s)
        for node1, node2 in self._graph.edges():
            print(f'"{node1.cell_id}" -- "{node2.cell_id}"', file=s)
        print("}", file=s)
        return s.getvalue()

    @property
    def graph(self):
        return self._graph

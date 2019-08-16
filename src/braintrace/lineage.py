"""
Lineage computation. A lineage is represented as a set of cells (CellSet).
"""
from typing import List, Dict
from collections import Counter, defaultdict
from io import StringIO
import math
import subprocess

from .graph import Graph
from .cell import Cell


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


# TODO this is unused
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

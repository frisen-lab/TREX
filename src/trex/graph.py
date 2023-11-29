from collections import OrderedDict


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

    def remove_node(self, node):
        for neighbor in list(self._nodes[node]):
            self._nodes[neighbor] = [n for n in self._nodes[neighbor] if n != node]
        del self._nodes[node]

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

    def induced_subgraph(self, nodes):
        nodes_set = set(nodes)
        new_nodes = {
            node: [neighbor for neighbor in self._nodes[node] if neighbor in nodes_set]
            for node in nodes
        }
        subgraph = Graph([])
        subgraph._nodes = new_nodes

        return subgraph

    def count_edges(self) -> int:
        """Return number of edges"""
        return sum(len(neighbors) for neighbors in self._nodes.values()) // 2

    def local_cut_vertices(self):
        """
        Return all vertices that, when removed, would lead to their neighborhood being split
        into two or more connected components.
        """
        vertices = []
        for node in self.nodes():
            neighbors = self.neighbors(node)
            if len(neighbors) < 2:
                continue
            subgraph = self.induced_subgraph(neighbors)
            if len(subgraph.connected_components()) > 1:
                vertices.append(node)

        return vertices

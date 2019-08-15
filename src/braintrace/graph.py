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

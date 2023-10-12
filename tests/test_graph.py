from trex.graph import Graph

import pytest

NODES = list("ABCDEFG")


@pytest.fixture
def graph():
    graph = Graph(NODES)
    #
    # (A) -- (B) -- (C)
    #  |    /  \
    #  |  /     \
    # (D)       (E)
    #
    # (F) -- (G)
    #
    for edge in ["AB", "BC", "AD", "DB", "BE", "FG"]:
        graph.add_edge(edge[0], edge[1])
    return graph


def test_nodes(graph):
    assert graph.nodes() == NODES


def test_edges(graph):
    edges = list(graph.edges())
    assert sorted(edges) == [("A", "B"), ("A", "D"), ("B", "C"), ("B", "D"), ("B", "E"), ("F", "G")]


def test_neighbors(graph):
    neighbors = graph.neighbors("B")
    assert sorted(neighbors) == ["A", "C", "D", "E"]


def test_remove_edge(graph):
    graph.remove_edge("C", "B")
    assert sorted(graph.neighbors("B")) == ["A", "D", "E"]


def test_connected_components(graph):
    components = graph.connected_components()
    assert len(components) == 2
    nodes1 = components[0].nodes()
    nodes2 = components[1].nodes()
    assert sorted(nodes1) == list("ABCDE")
    assert sorted(nodes2) == list("FG")


def test_remove_node(graph):
    graph.remove_node("B")
    assert sorted(graph.nodes()) == list("ACDEFG")
    assert sorted(graph.neighbors("A")) == ["D"]
    assert sorted(graph.edges()) == [("A", "D"), ("F", "G")]

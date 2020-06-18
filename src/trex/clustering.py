from typing import List, Callable, Set
from collections import defaultdict

from .graph import Graph


def kmers(s: str, k: int):
    """
    Yield all overlapping k-kmers of s.

    >>> list(kmers('hello', 3))
    ['hel', 'ell', 'llo']
    """
    for i in range(len(s) - k + 1):
        yield s[i:i+k]


def cluster_sequences(
        sequences: List[str], is_similar: Callable[[str, str], bool], k: int = 6
) -> List[List[str]]:
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

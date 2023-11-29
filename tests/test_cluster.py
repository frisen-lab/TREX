import textwrap

from trex.clustering import cluster_sequences
from trex.cli.run10x import is_similar


def parse_cluster(s):
    s = textwrap.dedent(s).strip()
    return [line.split()[0] for line in s.splitlines()]


def test_parse_cluster():
    assert (
        parse_cluster(
            """
        TGGGTCGTTATAA 1
        AGGGTCGTTATAA 1
        """
        )
        == ["TGGGTCGTTATAA", "AGGGTCGTTATAA"]
    )


def test_cluster_1():
    sequences = parse_cluster(
        """
        TGGGGGTCGAGTGATTAGAGGTCGTTATAA 1
        AGGGGGTCGAGTGAGTAGAGGTCGTTATAA 1
        AGGGGGTCGAGTGATTAGAGGTCGTTATAA 32
        AGGGGGTCGAGTGATTAGAGGTAGTTATAA 1
        -------CGAGTGATTAGAGGTCGTTATAA 2
        ---------AGTGATTAGAGGTCGTTATAA 1
        -----------------GAGGTCGTTATAA 1
        -------------------GGTCGTTATAA 1
        """
    )
    clusters = cluster_sequences(
        sequences,
        is_similar=lambda s, t: is_similar(s, t, min_overlap=10, max_hamming=5),
    )

    assert len(clusters) == 1
    assert sorted(clusters[0]) == sorted(sequences)

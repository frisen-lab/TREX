from trex.molecule import compute_consensus


def test_consensus_small():
    sequences = [
        "ACG",
        "ACG",
        "ATG",
    ]
    assert compute_consensus(sequences) == "ACG"
    assert compute_consensus(sequences[::-1]) == "ACG"


def test_consensus_with_gaps():
    sequences = [
        "ACG",
        "ACG",
        "-CG",
    ]
    assert compute_consensus(sequences) == "ACG"
    assert compute_consensus(sequences[::-1]) == "ACG"


def test_consensus_real_example():
    sequences = [
        "ACGGGGG-----------------------",
        "ACGGGGGTCAGGG-----------------",
        "ACGGGGGTCAGGGTAGCACGGCCGTTTCGC",
        "ACGGGGGTCAGGGTAGCACGGCCGTTTCGC",
        "ACGGGGGTCAGGGTAGCACGGCCGTTTCGC",
        "ACGGGGGTCAGGGTAGCACGGCCGTTTCGC",
        "ACGGGGGTCAGGGTAGCACGGCCGTTTCGC",
        "ACGGGGGTCAGGGTAGCACGGCCGTTTCGC",
        "ACGGGGGTCAGGGTAGCACGGCCGTTTCGC",
        "ACGGGGGTCAGGGTAGCACGGCCGTTTCGC",
        "ACGGGGGTCAGGGTAGCACGGCCGTTTCGC",
        "ACGGGGGTCAGGGTAGCACGGCCGTTTCGC",
        "---------------GCACGGCCGTTTCGC",
        "-----GGTCAGGGTAGCACGGCCGTTTCGC",
        "-----GGTCAGGGTAGCACGGCCGTTTCGC",
        "-----GGTCAGGGTAGCACGGCCGTTTCGC",
    ]
    assert compute_consensus(sequences) == "ACGGGGGTCAGGGTAGCACGGCCGTTTCGC"

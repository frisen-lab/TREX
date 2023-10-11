from trex.cli.run10x import is_similar

import pytest


def test_is_similar_unequal_length_raises():
    with pytest.raises(IndexError) as error:
        is_similar("ABC", "ABCDEF", 3, 1)
    error.match("same length")

    with pytest.raises(IndexError) as error:
        is_similar("A--BC", "A--BCDEF", 3, 1)
    error.match("same length")


@pytest.mark.parametrize("s,t,similar", [
    ("AAAAAAAAAA", "AAAAAAAAAA", True),
    ("CAAAAAAAAA", "AAAAAAAAAA", True),
    ("CCAAAAAAAA", "AAAAAAAAAA", False),

    ("", "", False),  # too little overlap
    ("ACGT", "ACGT", False),  # too little overlap
    ("TACGT", "TACGT", True),
    ("TACGT", "TATGT", True),
    ("TACGT", "TATGT", True),
    ("TACGT", "TTTGT", False),  # too many differences

    ("------ACGT", "------ACGT", False),  # too little overlap
    ("-----TACGT", "-----TACGT", True),
    ("-----TACGT", "-----TATGT", True),
    ("-----TACGT", "GGGGGTATGT", True),
    ("TACGT-----", "TATGT-----", True),
    ("TACGT-----", "TATGTGGGGG", True),
    ("TAC-----GT", "TAT-----GT", True),
    ("TAC-----GT", "TATGGGGGGT", True),
    ("--TAC---GT--", "--TATGGGGT--", True),
    ("--TAC---GTGG", "--TATGGGGT--", True),
    ("GGTAC---GT--", "GGTATGGGGT--", True),
])
def test_is_similar(s, t, similar):
    min_overlap = 5
    max_hamming = 1
    assert is_similar(s, t, min_overlap, max_hamming) == similar
    # Must be symmetric
    assert is_similar(t, s, min_overlap, max_hamming) == similar

from trex.cli.run10x import is_similar, SimilaritySet

import pytest


def test_is_similar_unequal_length_raises():
    with pytest.raises(IndexError) as error:
        is_similar("ABC", "ABCDEF", 3, 1)
    error.match("same length")

    with pytest.raises(IndexError) as error:
        is_similar("A--BC", "A--BCDEF", 3, 1)
    error.match("same length")


@pytest.mark.parametrize(
    "s,t,similar",
    [
        ("AAAAAAAAAA", "AAAAAAAAAA", True),
        ("CAAAAAAAAA", "AAAAAAAAAA", True),
        ("CCAAAAAAAA", "AAAAAAAAAA", False),
        ("", "", False),  # too little overlap
        ("ACGT", "ACGT", False),  # too little overlap
        ("TACGT", "TACGT", True),
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
        ("--TAC---TTGG", "--TATGGGGT--", False),
        ("GGTAC---GT--", "GGTATGGGGT--", True),
        ("---00TACGT", "GGGGGTATGT", True),  # Mixing 0 and -
        ("---00TACGT", "GGGGGTATGA", False),
    ],
)
def test_is_similar(s, t, similar):
    min_overlap = 5
    max_hamming = 1
    assert is_similar(s, t, min_overlap, max_hamming) == similar
    # Must be symmetric
    assert is_similar(t, s, min_overlap, max_hamming) == similar

    assert is_similar(s.replace("-", "0"), t, min_overlap, max_hamming) == similar
    assert is_similar(s, t.replace("-", "0"), min_overlap, max_hamming) == similar


@pytest.mark.parametrize(
    "s,strings,similar",
    [
        ("AAAAAAAAAA", ["AAAAAAAAAA", "TTTTTTTAAA"], True),
        ("TAAAAAAAAA", ["AAAAAAAAAA", "TTTTTTTAAA"], False),
        ("TTTTTTTAAA", ["AAAAAAAAAA", "TTTTTTTAAA"], True),
        ("TTTTTTTAAG", ["AAAAAAAAAA", "TTTTTTTAAA"], False),
        ("TTTTTT--AA", ["AAAAAAAAAA", "TTTTTTTAAA"], True),
        ("----AAAAAA", ["AAAAAAAAAA", "TTTTTTTAAA"], True),
        ("AAAAAAAAAA", ["AAAAAAAAAA"], True),
        ("AAAAAAAAAA", ["TTTTTTTAAA"], False),
    ],
)
def test_similarity_set(s, strings, similar):
    similarity_set = SimilaritySet(strings)
    assert similarity_set.contains(s) == similar

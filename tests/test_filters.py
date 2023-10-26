from trex.filters import is_low_complexity


def test_is_low_complexity():
    assert is_low_complexity("AAAAAAAAAAAAAAAAAAAAAAAAAAGAAA")
    assert is_low_complexity("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    assert is_low_complexity("AAAAA-------------------------")
    assert not is_low_complexity("ACGTACGTACGTACGTACGTACGTACGTAC")
    assert not is_low_complexity("ACGTACGTGGGGGGGGGGGGGGGGGGGGGGG")

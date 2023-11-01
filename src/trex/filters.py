from collections import Counter
from scipy.stats import entropy

# Threshold (in bits) for a low-complexity string (cloneID)
MIN_SHANNON_ENTROPY = 1.0


def is_low_complexity(s: str) -> bool:
    """
    Return True for low-complexity strings such as
    AAAAAAAAAAAAAAAAAAAAAAAAAAGAAA
    or even
    T-----------------------------
    """
    character_counts = list(Counter(s).values())
    return entropy(character_counts, base=2) < MIN_SHANNON_ENTROPY

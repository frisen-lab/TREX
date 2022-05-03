"""
A molecule is a DNA/RNA fragment that has potentially been sequenced multiple times
"""
from typing import NamedTuple, List
from collections import defaultdict
import numpy as np

from .bam import Read


class Molecule(NamedTuple):
    umi: str
    cell_id: str
    clone_id: str
    read_count: int


def compute_molecules(reads: List[Read]) -> List[Molecule]:
    """
    Group reads by cell ID and UMI into molecules.

    The cloneID of the molecule is the consensus of the cloneIDs in the group.
    """
    groups = defaultdict(list)
    for read in reads:
        groups[(read.umi, read.cell_id)].append(read.clone_id)

    molecules = []
    for (umi, cell_id), clone_ids in groups.items():
        clone_id_consensus = compute_consensus(clone_ids)
        molecules.append(
            Molecule(
                cell_id=cell_id,
                umi=umi,
                clone_id=clone_id_consensus,
                read_count=len(clone_ids),
            )
        )

    sorted_molecules = sorted(
        molecules, key=lambda mol: (mol.cell_id, mol.clone_id, mol.umi)
    )

    return sorted_molecules


def compute_consensus(sequences):
    """
    Compute a consensus for a set of sequences.

    All sequences must have the same length.
    """
    if len(sequences) == 1:
        return sequences[0]
    assert sequences

    # TODO
    # Ensure that the sequences are actually somewhat similar

    letters = np.array(["A", "C", "G", "T", "-", "0"])
    length = len(sequences[0])
    matrix = np.zeros([length, 6], dtype="float16")
    for sequence in sequences:
        align = np.zeros([length, 6], dtype="float16")
        for (i, ch) in enumerate(sequence):
            # turns each base into a number and position in numpy array
            if ch == "A":
                align[i, 0] = 1
            elif ch == "C":
                align[i, 1] = 1
            elif ch == "G":
                align[i, 2] = 1
            elif ch == "T":
                align[i, 3] = 1
            elif ch == "-":
                align[i, 4] = 0.1
            elif ch == "0":
                align[i, 5] = 0.1
        matrix += align

    # calculate base with maximum count for each position
    bin_consens = np.argmax(matrix, axis=1)
    # convert maximum counts into consensus sequence
    return "".join(letters[bin_consens])

"""
A molecule is a DNA/RNA fragment that has potentially been sequenced multiple times
"""
import re
from tinyalign import hamming_distance
from typing import NamedTuple, List
from collections import defaultdict, Counter
import numpy as np

from .bam import Read
from .clustering import cluster_sequences


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


def remove_odd_clone_ids(molecules: List[Molecule],
                        min_bases_detected: int = 7) -> List[Molecule]:
    """Discard CloneIDs with a less than min_bases_detected or a
    single (repeated) base detected."""
    def acceptable_clone_id(clone_id):
        detected_clone_id = re.sub("[-0]", "", clone_id)
        if len(detected_clone_id) <= min_bases_detected or len(set(detected_clone_id)) <= 1:
            return False
        else:
            return True

    new_molecules = [mol for mol in molecules if acceptable_clone_id(
        mol.clone_id)]
    return new_molecules


def correct_clone_ids_per_cell(
    molecules: List[Molecule], max_hamming: int, min_overlap: int = 20
) -> List[Molecule]:
    """
    Attempt to correct sequencing errors in the CloneID sequences of all
    molecules looking for similar sequences in the same cell.
    """
    # Obtain all cloneIDs (including those with '-' and '0')
    clone_ids = [m.clone_id for m in molecules]

    # Count the full-length cloneIDs
    counts = Counter(clone_ids)

    cell_dict = defaultdict(list)
    for molecule in molecules:
        cell_dict[molecule.cell_id].append((molecule.umi, molecule.clone_id))

    # Cluster them by Hamming distance
    def is_similar(s, t):
        # m = max_hamming
        bad_chars = {'-', '0'}
        if bad_chars & set(s) or bad_chars & set(t):
            # Remove suffix and/or prefix where sequences do not overlap
            s = s.lstrip("-0")
            t = t[-len(s):]
            t = t.lstrip("-0")
            s = s[-len(t):]
            s = s.rstrip("-0")
            if len(s) < min_overlap:
                return False
            t = t[: len(s)]
            t = t.rstrip("-0")
            if len(t) < min_overlap:
                return False
            s = s[: len(t)]
            # TODO allowed Hamming distance should be reduced relative to the
            #  overlap length
            # m = max_hamming * len(s) / len(original_length_of_s)
        return hamming_distance(s, t) <= max_hamming

    cell_correction_map = defaultdict(dict)

    for cell_id, cell_molecules in cell_dict.items():
        cell_clone_ids = set([m[1] for m in cell_molecules])

        if len(cell_clone_ids) > 1:
            cell_clone_ids = list(cell_clone_ids)
            clusters = cluster_sequences(cell_clone_ids, is_similar=is_similar,
                                         k=0)

            for cluster in clusters:
                if len(cluster) > 1:
                    # Pick most frequent cloneID as representative
                    longest = max([len(re.sub('[0-]', '', x)) for x in cluster])
                    subcluster = [x for x in cluster if
                                  len(re.sub('[0|-]', '', x)) == longest]
                    representative = max(subcluster,
                                         key=lambda bc: (counts[bc], bc))
                    cell_correction_map[cell_id].update(
                        {clone_id: representative for clone_id in cluster if
                         clone_id != representative})

    def get_correct_molecule(molecule):
        this_correction_map = cell_correction_map.get(molecule.cell_id, None)
        if this_correction_map is not None:
            new_clone_id = this_correction_map.get(molecule.clone_id,
                                                  molecule.clone_id)
            molecule = molecule._replace(clone_id=new_clone_id)
        return molecule

    # Create a new list of molecules in which the cloneIDs have been replaced
    # by their representatives
    return [get_correct_molecule(molecule) for molecule in molecules]

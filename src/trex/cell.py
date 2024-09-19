from dataclasses import dataclass
from typing import Dict, List
from collections import defaultdict, OrderedDict, Counter

from .molecule import Molecule


@dataclass(frozen=True, order=True)
class Cell:
    cell_id: str
    counts: Dict[str, int]  # UMI counts

    def __hash__(self):
        return hash(self.cell_id)


def compute_cells(
    sorted_molecules: List[Molecule], minimum_clone_id_length: int
) -> List[Cell]:
    """
    Group molecules by cell id
    """
    cell_id_groups = defaultdict(list)
    for molecule in sorted_molecules:
        clone_id = molecule.clone_id
        dashes = clone_id.count("-")
        zeroes = clone_id.count("0")
        if len(clone_id) - dashes - zeroes >= minimum_clone_id_length:
            cell_id_groups[molecule.cell_id].append(molecule)

    cells = []
    for cell_id, molecules in cell_id_groups.items():
        clone_ids = [molecule.clone_id for molecule in molecules]
        counts = OrderedDict(
            sorted(Counter(clone_ids).most_common(), key=lambda x: x[0].count("-"))
        )
        cells.append(Cell(cell_id=cell_id, counts=counts))
    return cells

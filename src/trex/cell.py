from typing import NamedTuple, Dict, List
from collections import defaultdict, OrderedDict, Counter

from .molecule import Molecule


class Cell(NamedTuple):
    cell_id: str
    counts: Dict[str, int]

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
        pure_li = clone_id.strip("-")
        # TODO may not work as intended (strip only removes prefixes and suffixes)
        pure_bc0 = clone_id.strip("0")
        if (
            len(pure_li) >= minimum_clone_id_length
            and len(pure_bc0) >= minimum_clone_id_length
        ):
            cell_id_groups[molecule.cell_id].append(molecule)

    cells = []
    for cell_id, molecules in cell_id_groups.items():
        clone_ids = [molecule.clone_id for molecule in molecules]
        counts = OrderedDict(
            sorted(Counter(clone_ids).most_common(), key=lambda x: x[0].count("-"))
        )
        cells.append(Cell(cell_id=cell_id, counts=counts))
    return cells

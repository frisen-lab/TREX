from typing import NamedTuple, Dict, List
from collections import defaultdict, OrderedDict, Counter

from .bamS3 import Read


class Cell(NamedTuple):
    cell_id: str
    read_counts: Dict[str, int]

    def __hash__(self):
        return hash(self.cell_id)


def compute_cells(sorted_reads: List[Read], minimum_clone_id_length: int) -> List[Cell]:
    """
    Group reads by cell id
    """
    cell_id_groups = defaultdict(list)
    for read in sorted_reads:
        clone_id = read.clone_id
        pure_li = clone_id.strip('-')
        # TODO may not work as intended (strip only removes prefixes and suffixes)
        pure_bc0 = clone_id.strip('0')
        if len(pure_li) >= minimum_clone_id_length and len(pure_bc0) >= minimum_clone_id_length:
            cell_id_groups[read.cell_id].append(read)

    cells = []
    for cell_id, reads in cell_id_groups.items():
        clone_ids = [read.clone_id for read in reads]
        read_counts = OrderedDict(sorted(Counter(clone_ids).most_common(),
            key=lambda x: x[0].count('-')))
        cells.append(Cell(cell_id=cell_id, read_counts=read_counts))
    return cells

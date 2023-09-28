import logging
from typing import List

import pandas as pd

from .cell import Cell
from .molecule import Molecule


class NiceFormatter(logging.Formatter):
    """
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).

    Based on http://stackoverflow.com/a/9218261/715090 .
    """

    def format(self, record):
        if record.levelno != logging.INFO:
            record.msg = "{}: {}".format(record.levelname, record.msg)
        return super().format(record)


def molecule_list_to_dataframe(molecules: List[Molecule]) -> pd.DataFrame:
    """Convert list of Molecules to pandas DataFrame"""
    molecule_dict = {"cell_id": [], "umi": [], "clone_id": [], "read_count": []}

    for molecule in molecules:
        molecule_dict["cell_id"].append(molecule.cell_id)
        molecule_dict["umi"].append(molecule.umi)
        molecule_dict["clone_id"].append(molecule.clone_id)
        molecule_dict["read_count"].append(getattr(molecule, "read_count", -1))

    return pd.DataFrame(molecule_dict)


def dataframe_to_molecule_list(df: pd.DataFrame) -> List[Molecule]:
    """Convert pandas DataFrame to list of Molecules"""
    molecules = [
        Molecule(
            umi=mol.umi,
            cell_id=mol.cell_id,
            clone_id=mol.clone_id,
            read_count=mol.get("read_count", -1),
        )
        for r, mol in df.iterrows()
    ]

    sorted_molecules = sorted(
        molecules, key=lambda mol: (mol.cell_id, mol.clone_id, mol.umi)
    )

    return sorted_molecules


def dataframe_to_cell_list(df: pd.DataFrame) -> List[Cell]:
    """Convert pandas DataFrame to list of Cells"""
    cell_list = []
    for r, cell in df.groupby("cell_id", observed=True):
        cell_id = r
        counts = {
            clone_id: counts
            for clone_id, counts in zip(cell.clone_id.values, cell.counts.values)
        }
        cell_list.append(Cell(cell_id=cell_id, counts=counts))

    return cell_list

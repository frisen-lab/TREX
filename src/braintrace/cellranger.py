"""
Dealing with CellRanger "outs/" files
"""
from pathlib import Path
from typing import Set

from xopen import xopen


class CellRangerError(Exception):
    pass


class CellRangerOuts2:
    """CellRanger 2 "outs/" directory structure"""
    MATRICES = 'filtered_gene_bc_matrices'
    BARCODES = 'barcodes.tsv'
    BAM = 'possorted_genome_bam.bam'

    def __init__(self, path: Path, genome_name: str = None):
        self.path = path
        matrices_path = path / self.MATRICES
        if not matrices_path.exists():
            raise CellRangerError(
                f"Directory '{self.MATRICES}/' must exist in the given outs/ directory")
        self.matrices_path: Path = matrices_path
        self.bam = path / self.BAM
        self.sample_dir = path.parent
        self.genome_dir = self._detect_genome_dir(genome_name)
        self.barcodes_path = self.genome_dir / self.BARCODES

    def _detect_genome_dir(self, genome_name):
        if genome_name is not None:
            return self.matrices_path / genome_name

        genomes = [p for p in self.matrices_path.iterdir() if p.is_dir()]
        if not genomes:
            raise CellRangerError(
                f"No subfolders found in the '{self.matrices_path}' folder")
        if len(genomes) > 1:
            message = "Exactly one genome folder expected in the " \
                f"'{self.matrices_path}' folder, but found:"
            for g in genomes:
                message += f'\n  {g!r}'
            raise CellRangerError(message)
        return genomes[0]

    def cellids(self) -> Set[str]:
        """
        Read barcodes.tsv, which contains a list of corrected and approved cellIDs like this:

        AAACCTGAGCGACGTA-1
        AAACCTGCATACTCTT-1
        """
        with xopen(self.barcodes_path) as f:
            ids = []
            for line in f:
                line = line.strip('\n')
                ids.append(line)
        return set(ids)


class CellRangerOuts3(CellRangerOuts2):
    MATRICES = 'filtered_feature_bc_matrix'
    BARCODES = 'barcodes.tsv.gz'

    def _detect_genome_dir(self, _genome_name):
        return self.matrices_path


def make_cellranger_outs(outs_path: Path, *args, **kwargs):
    """Detect CellRanger outs/ format and return an appropriate instance of CellRangerOuts2/3"""
    outs_path = Path(outs_path)
    if (outs_path / 'filtered_gene_bc_matrices').exists():
        return CellRangerOuts2(outs_path, *args, **kwargs)
    else:
        return CellRangerOuts3(outs_path, *args, **kwargs)

from pathlib import Path
import logging

from .bam import read_bam
from .cellranger import make_cellranger
from .error import BraintraceError


logger = logging.getLogger(__name__)


class DatasetReader:
    def __init__(self, output_dir: Path, genome_name, chromosome, start, end, prefix: bool):
        self.output_dir = output_dir
        self.genome_name = genome_name
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.prefix = prefix

    def read_one(self, path, output_bam_path):
        cellranger_dir = make_cellranger(path, self.genome_name)
        allowed_cell_ids = cellranger_dir.cellids()
        reads = read_bam(
            cellranger_dir.bam, output_bam_path, allowed_cell_ids, self.chromosome, self.start, self.end)
        assert len(reads) > 0, "No reads"
        return reads

    def read_all(self, transcriptome_inputs, amplicon_inputs, names, allowed_cell_ids):
        n_transcriptome = len(transcriptome_inputs)
        n_amplicon = len(amplicon_inputs)
        assert n_transcriptome > 0
        assert n_amplicon == 0 or n_amplicon == n_transcriptome
        assert n_transcriptome == len(names)

        # Cases:
        # - 1 transcriptome, 0 amplicon -- do not modify cell ids
        # - 1 transcriptome, 1 amplicon -- modify cell ids?
        # - n transcriptome, 0 amplicon -- modify cell ids!
        # - n transcriptome, n amplicon -- modify cell ids!

        if n_transcriptome == 1 and n_amplicon == 0:
            reads = self.read_one(transcriptome_inputs[0], self.output_dir / "entries.bam")
        elif n_transcriptome > 1 and n_amplicon == 0:
            datasets = []
            for path, name in zip(transcriptome_inputs, names):
                assert name is not None
                datasets.append(
                    self.read_one(path, self.output_dir / (name + "_entries.bam"))
                )
            reads = self.merge_datasets(datasets, names)
        else:
            assert n_amplicon > 0
            raise BraintraceError("Reading amplicon data not supported at the moment")

        if allowed_cell_ids:
            logger.debug("Allowed cell ids:\n- %s\n  ...", "\n- ".join(list(allowed_cell_ids)[:10]))
            logger.debug("Cell ids of reads:\n- %s\n  ...", "\n- ".join(r.cell_id for r in reads[:10]))
            reads = [r for r in reads if r.cell_id in allowed_cell_ids]
        if not reads:
            raise BraintraceError("No reads left after --filter-cellids filtering")
        return reads

    def merge_datasets(self, datasets, names):
        if self.prefix:
            def combine(cell_id, affix):
                return affix + "_" + cell_id
        else:
            def combine(cell_id, affix):
                return cell_id + "_" + affix
        reads = []
        for dataset, name in zip(datasets, names):
            for read in dataset:
                reads.append(read._replace(cell_id=combine(read.cell_id, name)))

        return reads

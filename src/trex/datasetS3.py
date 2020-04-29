from pathlib import Path
from itertools import zip_longest
import logging

from .bamS3 import read_bam



logger = logging.getLogger(__name__)


class DatasetReader:
    def __init__(self, output_dir: Path, chromosome, start, end, prefix: bool):
        self.output_dir = output_dir
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.prefix = prefix

    def read_one(self, path, output_bam_path):
        reads = read_bam(
            path,
            output_bam_path,
            self.chromosome,
            self.start,
            self.end,
        )
        assert len(reads) > 0, "No reads"
        return reads

    def read_all(self, transcriptome_inputs, amplicon_inputs, names, allowed_cell_ids):
        n_transcriptome = len(transcriptome_inputs)
        n_amplicon = len(amplicon_inputs)
        assert n_transcriptome > 0
        assert n_amplicon == 0 or n_amplicon == n_transcriptome
        assert n_transcriptome == len(names)

        if n_transcriptome == 1:
            reads = self.read_one(transcriptome_inputs[0], self.output_dir / "entries.bam")
            if n_amplicon == 1:
                reads.extend(
                    self.read_one(amplicon_inputs[0], self.output_dir / "amplicon_entries.bam")
                )
        else:
            datasets = []
            for *paths, name in zip_longest(transcriptome_inputs, amplicon_inputs, names):
                assert name is not None
                reads = self.read_one(paths[0], self.output_dir / (name + "_entries.bam"))
                if paths[1]:
                    reads.extend(
                        self.read_one(paths[1], self.output_dir / (name + "_amplicon_entries.bam"))
                    )
                datasets.append(reads)
            reads = self.merge_datasets(datasets, names)

        if allowed_cell_ids:
            logger.debug("Allowed cell ids:\n- %s\n  ...", "\n- ".join(list(allowed_cell_ids)[:10]))
            logger.debug(
                "Cell ids of reads:\n- %s\n  ...", "\n- ".join(r.cell_id for r in reads[:10]))
            reads = [r for r in reads if r.cell_id in allowed_cell_ids]
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

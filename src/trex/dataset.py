from pathlib import Path
from itertools import zip_longest
import logging
from typing import Optional, List

from .bam import read_bam, write_outbam
from .cellranger import make_cellranger


logger = logging.getLogger(__name__)


def find_files(path, cell_id_tag):
    if cell_id_tag == "CB":
        files = [path]
    else:
        files = sorted(Path(path).glob("*.bam"))
        if len(files) == 0:
            files = [path]
    return files


class DatasetReader:
    def __init__(
        self, output_dir: Path, genome_name, chromosome, start, end, prefix: bool
    ):
        self.output_dir = output_dir
        self.genome_name = genome_name
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.prefix = prefix

    def read_multiple(self, input, output_bam_path, cell_id_tag="CB", require_umis=True):
        allowed_cell_ids = None
        files = find_files(input, cell_id_tag)
        all_reads = []
        all_reads_seq = []
        for file in files:
            bam = file
            if cell_id_tag == "CB":
                cellranger_dir = make_cellranger(file, self.genome_name)
                allowed_cell_ids = cellranger_dir.cellids()
                bam = cellranger_dir.bam

            reads, reads_seq, input_bam_path = read_bam(
                bam,
                allowed_cell_ids,
                self.chromosome,
                self.start,
                self.end,
                require_umis,
                cell_id_tag,
            )
            all_reads.extend(reads)
            all_reads_seq.extend(reads_seq)

        write_outbam(all_reads_seq=all_reads_seq, output_bam_path=output_bam_path,
                     input_bam_path=input_bam_path)

        return all_reads

    def read_all(
        self,
        transcriptome_inputs,
        amplicon_inputs,
        names: Optional[List[str]],
        allowed_cell_ids: Optional[List[str]],
        require_umis=True,
        cell_id_tag="CB",
    ):
        n_transcriptome = len(transcriptome_inputs)
        n_amplicon = len(amplicon_inputs)
        if names is None:
            names = [None] * n_transcriptome
        assert n_transcriptome > 0
        assert n_amplicon == 0 or n_amplicon == n_transcriptome
        assert n_transcriptome == len(names)

        if n_transcriptome == 1:
            reads = self.read_multiple(
                transcriptome_inputs[0],
                self.output_dir / "entries.bam",
                cell_id_tag,
                require_umis,
            )
            if n_amplicon == 1:
                reads.extend(
                    self.read_multiple(
                        amplicon_inputs[0],
                        self.output_dir / "amplicon_entries.bam",
                        cell_id_tag,
                        require_umis,
                    )
                )
        else:
            datasets = []
            for *paths, name in zip_longest(
                transcriptome_inputs, amplicon_inputs, names
            ):
                assert name is not None
                reads = []
                reads.exend(
                    self.read_multiple(
                        paths[0],
                        self.output_dir / (name + "_entries.bam"),
                        cell_id_tag,
                        require_umis,
                    )
                )
                if paths[1]:
                    reads.extend(
                            self.read_multiple(
                                paths[1],
                                self.output_dir / (name + "_amplicon_entries.bam"),
                                cell_id_tag,
                                require_umis,
                            )
                        )
                datasets.append(reads)
            reads = self.merge_datasets(datasets, names)

        if allowed_cell_ids:
            logger.debug(
                "Allowed cell ids:\n- %s\n  ...",
                "\n- ".join(list(allowed_cell_ids)[:10]),
            )
            logger.debug(
                "Cell ids of reads:\n- %s\n  ...",
                "\n- ".join(r.cell_id for r in reads[:10]),
            )
            reads = [r for r in reads if r.cell_id in allowed_cell_ids]

        if require_umis:
            umis = list()
            for read in reads:
                if read.umi:
                    umis.append(read.umi)
            assert len(umis) > 0, "No UMIs"
        assert len(reads) > 0, "No reads"

        sorted_reads = sorted(reads, key=lambda rd: (rd.cell_id, rd.clone_id))
        return sorted_reads

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

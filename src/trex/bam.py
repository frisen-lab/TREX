"""
Extract reads from BAM files
"""
import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple
from typing import Optional

from pysam import AlignmentFile, AlignedSegment
import pysam

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class Read:
    umi: Optional[str]
    cell_id: str
    clone_id: str


def read_bam(
    bam_path: Path,
    allowed_cell_ids: List[str],
    chr_name: str,
    clone_id_start=None,
    clone_id_end=None,
    require_umis=True,
    cell_id_tag="CB",
):
    """
    bam_path -- path to input BAM file or folder with input bam files
    output_dir -- path to an output directory into which a BAM file is written that contais all
        reads on the chromosome that have the required tags.
    """
    with AlignmentFile(bam_path) as alignment_file:
        if chr_name is None:
            chr_name = alignment_file.references[-1]

        if clone_id_start is None or clone_id_end is None:
            if clone_id_start is not None or clone_id_end is not None:
                raise ValueError(
                    "Either both or none of cloneID start and end must be provided"
                )
            clone_id_start, clone_id_end = detect_clone_id_location(
                alignment_file, chr_name
            )
            if clone_id_end - clone_id_start < 10:
                raise ValueError("Auto-detected cloneID too short, something is wrong")
        logger.info(
            f"Reading cloneIDs from {chr_name}:{clone_id_start + 1}-{clone_id_end} "
            f"in {bam_path}"
        )
        reads, no_umi, no_cell_id, reads_seq, alignment_file = read_alignment_file(
            alignment_file,
            allowed_cell_ids,
            chr_name,
            clone_id_start,
            clone_id_end,
            require_umis,
            cell_id_tag,
        )
    if require_umis:
        logger.info(
            f"Found {len(reads)} reads with usable cloneIDs and UMIs. "
            f"Skipped {no_cell_id} without cell id, "
            f"{no_umi} without UMI."
        )
    else:
        logger.info(
            f"Found {len(reads)} reads with usable cloneIDs. Skipped {no_cell_id} without cell id"
        )

    return reads, reads_seq, bam_path


def read_alignment_file(
    alignment_file: AlignmentFile,
    allowed_cell_ids: List[str],
    chr_name: str,
    clone_id_start: int,
    clone_id_end: int,
    require_umis=True,
    cell_id_tag="CB",
) -> Tuple[List[Read], int, int]:
    """
    bam_path -- path to input BAM file
    output_dir -- path to an output directory into which a BAM file is written that contais all
        reads on the chromosome that have the required tags.
    """

    clone_id_extractor = CachedCloneIdExtractor(clone_id_start, clone_id_end)
    # Fetches those reads aligning to the artifical, clone-id-containing chromosome
    reads = []
    reads_seq = []
    no_cell_id = no_umi = 0
    start, stop = max(0, clone_id_start - 10), clone_id_end + 10
    for read in alignment_file.fetch(chr_name, start, stop):
        # Collect all read sequences for output BAM file
        reads_seq.append(read)
        # Skip reads without cellID or UMI
        has_cell_id = read.has_tag(cell_id_tag)
        has_umi = read.has_tag("UB")
        if has_umi:
            if len(read.get_tag("UB")) == 0:
                no_umi += 1
        if not has_cell_id:
            no_cell_id += 1
        if not has_umi:
            no_umi += 1
        if not has_cell_id or (require_umis and not has_umi):
            continue
        cell_id = read.get_tag(cell_id_tag)
        if not isinstance(cell_id, str):
            raise ValueError(f"{cell_id_tag} tag must be a string")
        if allowed_cell_ids and cell_id not in allowed_cell_ids:
            no_cell_id += 1
            continue

        umi: Optional[str] = read.get_tag("UB")
        if umi == "":
            umi = None
        if cell_id_tag == "CB":
            if not cell_id.endswith("-1"):
                raise ValueError(
                    f"A cell id ({cell_id!r}) was found that does not end in '-1'. "
                    "Currently, this type of data cannot be used"
                )
            cell_id = cell_id[:-2]
        clone_id = clone_id_extractor.extract(read)
        if clone_id is None:
            # Read does not cover the cloneID
            continue
        if require_umis and not umi:
            # Read does not have a UMI
            continue
        reads.append(Read(cell_id=cell_id, umi=umi, clone_id=clone_id))

    logger.debug(
        f"CloneID extractor cache hits: {clone_id_extractor.hits}. "
        f"Cache misses: {clone_id_extractor.misses}"
    )
    return reads, no_umi, no_cell_id, reads_seq, alignment_file


def write_outbam(
    all_reads_seq: AlignedSegment, output_bam_path: str, input_bam_path: str
):
    # Write the passing alignments to a separate file
    print("I am writing an output file")

    alignment_file = pysam.AlignmentFile(input_bam_path, "rb")

    with AlignmentFile(output_bam_path, "wb", template=alignment_file) as out_bam:
        for read in all_reads_seq:
            out_bam.write(read)


class CachedCloneIdExtractor:
    """
    Caching cloneID extraction results helps when there are many identical reads
    (amplicon data)
    """

    def __init__(self, clone_id_start, clone_id_end):
        self._cache = dict()
        self._start = clone_id_start
        self._end = clone_id_end
        self.hits = 0
        self.misses = 0
        self._prev_start = None

    def extract(self, read):
        if self._prev_start != read.reference_start:
            self._cache = dict()
            self._prev_start = read.reference_start
        cache_key = (read.reference_start, read.query_sequence)
        try:
            clone_id = self._cache[cache_key]
            self.hits += 1
        except KeyError:
            self.misses += 1
            clone_id = self._extract(read)
            self._cache[cache_key] = clone_id
        return clone_id

    def _extract(self, read):
        query_align_end = read.query_alignment_end
        query_align_start = read.query_alignment_start
        query_sequence = read.query_sequence
        # Extract cloneID
        clone_id = ["-"] * (self._end - self._start)
        bases = 0
        for query_pos, ref_pos in read.get_aligned_pairs():
            # Replace soft-clipping with an ungapped alignment extending into the
            # soft-clipped region, assuming the clipping occurred because the cloneID
            # region was encountered
            if ref_pos is None:
                # Soft clip or insertion
                if query_align_end <= query_pos:
                    # We are in a soft-clipped region at the 3' end of the read
                    ref_pos = read.reference_end + (query_pos - query_align_end)
                elif query_align_start > query_pos:
                    # We are in a soft-clipped region at the 5' end of the read
                    ref_pos = read.reference_start - (query_align_start - query_pos)
                # ref_pos remains None if this is an insertion

            if ref_pos is not None and self._start <= ref_pos < self._end:
                if query_pos is None:
                    # Deletion or intron skip
                    query_base = "0"
                else:
                    # Match or mismatch
                    query_base = query_sequence[query_pos]
                    bases += 1
                clone_id[ref_pos - self._start] = query_base

        if bases == 0:
            return None
        else:
            return "".join(clone_id)


def detect_clone_id_location(
    alignment_file: AlignmentFile, reference_name: str
) -> Tuple[int, int]:
    """
    Detect where the cloneID is located on the reference by inspecting the alignments.

    Return (clone_id_start, clone_id_end)
    """
    # Look for reference positions at which reads are soft-clipped at their 3' end
    starts = Counter()
    reference_length = alignment_file.get_reference_length(reference_name)
    for alignment in alignment_file.fetch(reference_name):
        clip_right = alignment.query_length - alignment.query_alignment_end
        if clip_right >= 5:
            starts[alignment.reference_end] += 1

    for clone_id_start, freq in starts.most_common(5):
        # Soft-clipping at the 5' end cannot be used to find the cloneID end when
        # the cloneID region is too far at the 3' end of the contig. Instead,
        # look at pileups and check base frequencies (over the cloneID, bases should
        # be roughly uniformly distributed).
        if clone_id_start >= reference_length:
            # The most common reference position that is soft clipped is often the 3' end
            # of the contig. Skip that.
            continue
        clone_id_end = clone_id_start
        for column in alignment_file.pileup(reference_name, start=clone_id_start):
            if column.reference_pos < clone_id_start:
                # See pileup() documentation
                continue
            bases = [
                p.alignment.query_sequence[p.query_position]
                for p in column.pileups
                if p.query_position is not None
            ]
            counter = Counter(bases)
            # Check whether one base dominates
            if counter.most_common()[0][1] / len(bases) > 0.95:
                # We appear to have found the end of the cloneID
                clone_id_end = column.reference_pos
                break
        if clone_id_end - clone_id_start >= 5:
            # Good enough
            return (clone_id_start, clone_id_end)
    raise ValueError(
        f"Could not detect cloneID location on chromosome {reference_name}"
    )

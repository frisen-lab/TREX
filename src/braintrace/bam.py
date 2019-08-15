"""
Extract reads from BAM files
"""
import logging
from collections import Counter
from pathlib import Path
from typing import NamedTuple

from pysam import AlignmentFile

logger = logging.getLogger(__name__)


class Read(NamedTuple):
    umi: str
    cell_id: str
    clone_id: str


def read_bam(bam_path: Path, output_dir: Path, allowed_cell_ids, chr_name, clone_id_start=None, clone_id_end=None, file_name_suffix="_entries", cellid_suffix=None):
    """
    bam_path -- path to input BAM file
    output_dir -- path to an output directory into which a BAM file is written that contais all
        reads on the chromosome that have the required tags.
    """
    cache_hit = cache_miss = 0
    with AlignmentFile(bam_path) as alignment_file:
        if chr_name is None:
            chr_name = alignment_file.references[-1]

        # TODO move out of here
        if clone_id_start is None or clone_id_end is None:
            if clone_id_start is not None or clone_id_end is not None:
                raise ValueError('Either both or none of clone id start and end must be provided')
            clone_id_start, clone_id_end = detect_clone_id_location(alignment_file, chr_name)
        logger.info(f"Reading clone ids from {chr_name}:{clone_id_start + 1}-{clone_id_end} "
            f"in {bam_path}")
        if clone_id_end - clone_id_start < 10:
            raise ValueError('Auto-detected clone id too short, something is wrong')
        output_bam_path = output_dir / (chr_name + file_name_suffix + '.bam')

        with AlignmentFile(output_bam_path, 'wb', template=alignment_file) as out_bam:
            # Fetches those reads aligning to the artifical, clone-id-containing chromosome
            reads = []
            unknown_ids = no_cell_id = no_umi = 0
            prev_start = None
            for read in alignment_file.fetch(chr_name, max(0, clone_id_start - 10), clone_id_end + 10):
                # Skip reads without cellID or UMI
                if not read.has_tag('CB') or not read.has_tag('UB'):
                    if not read.has_tag('CB'):
                        no_cell_id += 1
                    if not read.has_tag('UB'):
                        no_umi += 1
                    continue
                # Filters out reads that have not approved cellIDs
                cell_id = read.get_tag('CB')
                if cell_id not in allowed_cell_ids:
                    unknown_ids += 1
                    continue
                # Gives cellIDs an experiment-specific suffix
                if cellid_suffix is not None:
                    if cell_id.endswith("-1"):
                        cell_id = cell_id[:-2]
                    cell_id += cellid_suffix

                # Use a cache of clone id extraction results. This helps when there are
                # many identical reads (amplicon data)
                if prev_start != read.reference_start:
                    clone_id_cache = dict()
                    prev_start = read.reference_start
                cache_key = (read.reference_start, read.query_sequence)
                try:
                    clone_id = clone_id_cache[cache_key]
                    cache_hit += 1
                except KeyError:
                    cache_miss += 1
                    clone_id = extract_clone_id(clone_id_start, clone_id_end, read)
                    clone_id_cache[cache_key] = clone_id
                if clone_id is None:
                    # Read does not cover the clone id
                    continue
                reads.append(Read(cell_id=cell_id, umi=read.get_tag('UB'), clone_id=clone_id))

                # Write the passing alignments to a separate file
                out_bam.write(read)

            logger.info(f'Skipped {unknown_ids} reads with unrecognized cell ids '
                        f'(and {no_umi+no_cell_id} without UMI or cell id)')
    logger.info(f"Cache hits: {cache_hit}. Cache misses: {cache_miss}")
    sorted_reads = sorted(reads, key=lambda read: (read.umi, read.cell_id, read.clone_id))
    assert len(sorted_reads) == 0 or len(sorted_reads[0].clone_id) == clone_id_end - clone_id_start
    return sorted_reads


def detect_clone_id_location(alignment_file, reference_name):
    """
    Detect where the clone id is located on the reference by inspecting the alignments.

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
        # Soft-clipping at the 5' end cannot be used to find the clone id end when
        # the clone id region is too far at the 3' end of the contig. Instead,
        # look at pileups and check base frequencies (over the clone id, bases should
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
            bases = [p.alignment.query_sequence[p.query_position] for p in column.pileups if p.query_position is not None]
            counter = Counter(bases)
            # Check whether one base dominates
            if counter.most_common()[0][1] / len(bases) > 0.95:
                # We appear to have found the end of the clone id
                clone_id_end = column.reference_pos
                break
        if clone_id_end - clone_id_start >= 5:
            # Good enough
            return (clone_id_start, clone_id_end)
    raise ValueError(f'Could not detect clone id location on chromosome {reference_name}')


def extract_clone_id(clone_id_start, clone_id_end, read):
    query_align_end = read.query_alignment_end
    query_align_start = read.query_alignment_start
    query_sequence = read.query_sequence
    # Extract clone id
    clone_id = ['-'] * (clone_id_end - clone_id_start)
    bases = 0
    for query_pos, ref_pos in read.get_aligned_pairs():
        # Replace soft-clipping with an ungapped alignment extending into the
        # soft-clipped region, assuming the clipping occurred because the clone id
        # region was encountered.
        if ref_pos is None:
            # Soft clip or insertion
            if query_align_end <= query_pos:
                # We are in a soft-clipped region at the 3' end of the read
                ref_pos = read.reference_end + (query_pos - query_align_end)
            elif query_align_start > query_pos:
                # We are in a soft-clipped region at the 5' end of the read
                ref_pos = read.reference_start - (query_align_start - query_pos)
            # ref_pos remains None if this is an insertion

        if ref_pos is not None and clone_id_start <= ref_pos < clone_id_end:
            if query_pos is None:
                # Deletion or intron skip
                query_base = '0'
            else:
                # Match or mismatch
                query_base = query_sequence[query_pos]
                bases += 1
            clone_id[ref_pos - clone_id_start] = query_base

    if bases == 0:
        return None
    else:
        return "".join(clone_id)

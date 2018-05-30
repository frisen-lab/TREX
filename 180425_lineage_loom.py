"""
Description: Program for the extraction and filtering of random barcodes from single-cell
sequencing data

Preparation: Program processes cell ranger output files. Run cell ranger before. Follow
instructions => cellranger_instructions.sh

Run: Run program in cellranger 'outs' directory OR indicate path to 'outs'-directory via --path flag
"""

import os.path
import argparse
from collections import Counter
from collections import defaultdict
from collections import namedtuple
import operator
import warnings

import numpy as np
import pysam
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', 'Conversion of the second argument of issubdtype')
    import loompy


__author__ = 'leonie.von.berlin@stud.ki.se'


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome-name',
        help='Name of the genome as indicated in cell ranger count run with the flag --genome. '
             'Default: %(default)s',
        default='hg38_Tomato-N')
    parser.add_argument('--chromosome', '--chr',
        help='Barcode chromosome name as indicated in .fasta file. See cellranger_instructions.sh'
             'Default: %(default)s.',
        default='chrTomato-N')
    parser.add_argument('--path', '-p',
        help='Path to cell ranger "outs" directory. Default: current directory',
        default=os.getcwd())
    parser.add_argument('--name', '-n',
        help='name of the run and directory created by program. Default: lineage_run',
        default='lineage_run')
    parser.add_argument('--start', '-s',
        help='Position of first base INSIDE the barcode (with first base on position 0). '
             'Default: %(default)s',
        type=int, default=694)
    parser.add_argument('--end', '-e',
        help='Position of last base INSIDE the barcode (with first base on position 0). '
             'Default: %(default)s',
        type=int, default=724)
    parser.add_argument('--min-length', '-m',
        help='Minimum number of bases a barcode must have. Default: 10', type=int, default=10)
    parser.add_argument('--hamming',
        help='Minimum hamming distance allowed for two barcodes to be called similar. '
             'Default: %(default)s',
        type=int, default=4)
    parser.add_argument('-l', '--loom',
        help='If given, create loom-file from cell ranger and barcode data. '
             'File will have the same name as the run',
        action='store_true')
    return parser.parse_args()


def read_cellid_barcodes(path):
    """
    Read barcodes.tsv, which contains a list of corrected and approved cellIDs like this:

    AAACCTGAGCGACGTA-1
    AAACCTGCATACTCTT-1
    """
    with open(path) as f:
        ids = []
        for line in f:
            line = line.strip('\n')
            ids.append(line)
    return set(ids)


Read = namedtuple('Read', ['cell_id', 'umi', 'barcode'])

Molecule = namedtuple('Molecule', ['cell_id', 'umi', 'barcode', 'read_count'])

Cell = namedtuple('Cell', ['cell_id', 'barcode_counts'])


def read_bam(bam_path, output_bam_path, cell_ids, chr_name, barcode_start, barcode_end):
    """
    bam_path -- path to input BAM file
    output_bam_path -- path to an output BAM file. All reads on the chromosome that have the
        required tags are written to this file
    """
    alignment_file = pysam.AlignmentFile(bam_path)
    out_bam = pysam.AlignmentFile(output_bam_path, 'wb', template=alignment_file)

    # Fetches those reads aligning to the artifical, barcode-containing chromosome
    reads = []
    for read in alignment_file.fetch(chr_name):
        # Skip reads without cellID or UMI
        if not read.has_tag('CB') or not read.has_tag('UB'):
            continue
        # Filters out reads that have not approved cellIDs
        cell_id = read.get_tag('CB')
        if cell_id not in cell_ids:
            continue

        # Write the passing alignments to a separate file
        # TODO write only the alignments that actually cover the barcode region
        out_bam.write(read)

        query_align_end = read.query_alignment_end
        query_align_start = read.query_alignment_start

        # Extract barcode
        barcode = ['-'] * (barcode_end - barcode_start)
        for query_pos, ref_pos in read.get_aligned_pairs():
            # We replace soft-clipping with an ungapped alignment extending into the soft-clipped
            # regions, assuming the clipping occurred because the barcode region was encountered.
            if ref_pos is None:
                # Soft clip or insertion
                if query_align_end <= query_pos:
                    # We are in a soft-clipped region at the 3' end of the read
                    ref_pos = read.reference_end + (query_pos - query_align_end)
                elif query_align_start > query_pos:
                    # We are in a soft-clipped region at the 5' end of the read
                    ref_pos = read.reference_start - (query_align_start - query_pos)
                # ref_pos remains None if this is an insertion

            if ref_pos is not None and barcode_start <= ref_pos < barcode_end:
                if query_pos is None:
                    # Deletion or intron skip
                    query_base = '0'
                else:
                    # Match or mismatch
                    query_base = read.query_sequence[query_pos]
                barcode[ref_pos - barcode_start] = query_base

        barcode = ''.join(barcode)
        reads.append(Read(cell_id=cell_id, umi=read.get_tag('UB'), barcode=barcode))

    sorted_reads = sorted(reads, key=lambda read: (read.umi, read.cell_id, read.barcode))
    alignment_file.close()
    out_bam.close()
    assert len(sorted_reads) == 0 or len(sorted_reads[0].barcode) == barcode_end - barcode_start
    return sorted_reads


def compute_consensus(sequences):
    """
    Compute a consensus for a set of sequences.

    All sequences must have the same length.
    """
    if len(sequences) == 1:
        return sequences[0]
    assert sequences

    letters = np.array(['A', 'C', 'G', 'T', '-', '0'])
    length = len(sequences[0])
    consens_np = np.zeros([length, 6], dtype='float16')
    for sequence in sequences:
        align = np.zeros([length, 6], dtype='float16')
        for (i, ch) in enumerate(sequence):
            # turns each base into a number and position in numpy array
            if ch == 'A':
                align[i, 0] = 1
            elif ch == 'C':
                align[i, 1] = 1
            elif ch == 'G':
                align[i, 2] = 1
            elif ch == 'T':
                align[i, 3] = 1
            elif ch == '-':
                align[i, 4] = 0.1
            elif ch == '0':
                align[i, 5] = 0.1
        consens_np += align
    # calculate base with maximum count for each position
    bin_consens = np.argmax(align, axis=1)
    # convert maximum counts into consensus sequence
    return ''.join(letters[bin_consens])


def compute_molecules(sorted_reads):
    """
    - Forms groups of reads with identical CellIDs and UMIs => belong to one molecule

    - forms consensus sequence of all barcodes of one group,
    """
    groups = defaultdict(list)
    for read in sorted_reads:
        groups[(read.umi, read.cell_id)].append(read.barcode)

    molecules = []
    for (umi, cell_id), barcodes in groups.items():
        barcode_consensus = compute_consensus(barcodes)
        molecules.append(
            Molecule(cell_id=cell_id, umi=umi, barcode=barcode_consensus, read_count=len(barcodes)))

    sorted_molecules = sorted(molecules, key=lambda mol: (mol.cell_id, mol.barcode, mol.umi))

    return sorted_molecules


def compute_cells(sorted_molecules, minimum_barcode_length, minham):
    # 1. Forms groups of molecules (with set barcode minimum length) that have identical cellIDs
    #    => belong to one cell,
    # 2. counts number of appearances of each barcode in each group,
    # 3. starting from the barcode with the lowest count, compares to barcodes starting with
    #    the highest counts of a group and calculates hamming distance. If distance is below
    #    threshold, the two barcodes and counts are merged. Repetition until no barcodes with
    #    hamming distance below threshold can be found (note that this way of merging is greedy),
    # 4. Outputs for each cells all its barcodes and corresponding counts

    cell_id_groups = defaultdict(list)
    for molecule in sorted_molecules:
        barcode = molecule.barcode
        pure_bc = barcode.strip('-')
        # TODO may not work as intended (strip only removes prefixes and suffixes)
        pure_bc0 = barcode.strip('0')
        if len(pure_bc) >= minimum_barcode_length and len(pure_bc0) >= minimum_barcode_length:
            cell_id_groups[molecule.cell_id].append(molecule)

    cells = []

    # merges barcodes and counts below hamming distance
    found = False
    for cell_id, molecules in cell_id_groups.items():
        barcodes = [molecule.barcode for molecule in molecules]
        results = defaultdict(int)
        # barcodes sorted by counts
        mc = sorted(Counter(barcodes).most_common(), key=lambda x: len(x[0].strip('-')),
            reverse=True)
        while True:
            barcode, n = mc.pop(-1)  # takes out and remove barcode with lowest count from list
            if len(mc) == 0:  # or '0' in x: #if barcode is the last in the list or it contains insertions/deletions (cannot be compared) just keeps barcode without merging
                results[barcode] += n
                break
            for i, m in mc:  # goes through remaining barcodes in list
                hamming = 0
                overlap_count = 0
                for l, k in zip(barcode, i):
                    if l != '-' and k != '-':  # only compares base-containing and not empty position
                        overlap_count += 1  # counts the overlap of two barcodes
                        if l != k:
                            hamming += 1  # calculates hamming distance based on the similarity of each base-pair
                if hamming < minham and overlap_count != 0:  # filters out barcode-pairs with hamming distance below set threshold or with no overlapp
                    if i.count('-') == 0:  # only full barcodes are merged with other groups
                        results[i] += n
                    else:
                        results[barcode] += n
                    found = True
                    break

            if not found:  # barcodes that never undergo the hamming distance threshold, are not merged
                results[barcode] += n
            else:
                found = False

        cells.append(Cell(cell_id=cell_id, barcode_counts=results))
    return cells


def filter_cells(cells, molecules):
    """
    Filter barcodes according to two criteria:

    - Barcodes that have only a count of one and can be found in another cell are most
      likely results of contamination and are removed,
    - Barcodes that have only a count of one and are also only based on one read are
      also removed
    """
    overall_barcode_counts = Counter()
    for cell in cells:
        overall_barcode_counts.update(cell.barcode_counts)

    single_barcode_cell_ids = set()
    for molecule in molecules:
        if molecule.read_count == 1:
            single_barcode_cell_ids.add(molecule.cell_id)
    # or:
    # single_barcode_cell_ids = {m.cell_id for m in molecules if m.read_count == 1}

    # filters out barcodes with a count of one that appear in another cell
    new_cells = []
    for cell in cells:
        barcode_counts = cell.barcode_counts.copy()
        for barcode, count in cell.barcode_counts.items():
            if count > 1:
                # This barcode occurs more than once in this cell - keep it
                continue
            if overall_barcode_counts[barcode] > 1:
                # This barcode occurs also in other cells - remove it
                del barcode_counts[barcode]
            elif cell.cell_id in single_barcode_cell_ids:
                # FIXME
                # This is probably incorrect (too broad): The barcode is discarded if *any*
                # barcode in this cell is based on a single read
                del barcode_counts[barcode]
        new_cells.append(Cell(cell_id=cell.cell_id, barcode_counts=barcode_counts))
    return new_cells


def write_loom(cell_col, input_dir, run_name, len_bc):
    bc_dict = {'1': [], '2': [], '3': [], '4': [], '5': [], '6': []}
    cnt_dict = {'1': [], '2': [], '3': [], '4': [], '5': [], '6': []}
    cellid1 = []

    # brings the barcode data into a format where the most abundant barcode of the cells are in
    # one list, the second most abundant in another and so on. The same with counts
    for i in range(0, len(cell_col), 2):
        sort_d = sorted(cell_col[i + 1].items(), key=operator.itemgetter(1))
        sort_d.reverse()
        if len(sort_d) != 0:
            cellid1.append(cell_col[i])
            for j in range(0, 6):
                k = j + 1
                if j <= len(sort_d) - 1:
                    bc_dict[str(k)].append(sort_d[j][0])
                    cnt_dict[str(k)].append(sort_d[j][1])
                else:
                    bc_dict[str(k)].append('-')
                    cnt_dict[str(k)].append(0)

    # creates the loom file based on cellranger output files

    loom_name = os.path.basename(input_dir[:-5])
    loom_path = os.path.join(run_name, loom_name + '.loom')
    if not os.path.exists(loom_path):
        loompy.create_from_cellranger(input_dir[:-5], run_name)
    # connects to the just created loom file in order to modify it
    ds = loompy.connect(loom_path)
    # gets a list of all cellIDs appearing in the loom file
    all_cellIDs = ds.ca.CellID

    # brings barcode data into correct format for loom file.
    # Array must have same shape as all_cellIDs
    bc_fulldict = {'1': [], '2': [], '3': [], '4': [], '5': [], '6': []}
    cnt_fulldict = {'1': [], '2': [], '3': [], '4': [], '5': [], '6': []}
    for id1 in all_cellIDs:
        found = False
        for id2 in cellid1:
            if id1[(len(loom_name) + 1):] == id2:
                found = True
                index = cellid1.index(id2)
                bc_fulldict['1'].append(bc_dict['1'][index])
                bc_fulldict['2'].append(bc_dict['2'][index])
                bc_fulldict['3'].append(bc_dict['3'][index])
                bc_fulldict['4'].append(bc_dict['4'][index])
                bc_fulldict['5'].append(bc_dict['5'][index])
                bc_fulldict['6'].append(bc_dict['6'][index])

                cnt_fulldict['1'].append(cnt_dict['1'][index])
                cnt_fulldict['2'].append(cnt_dict['2'][index])
                cnt_fulldict['3'].append(cnt_dict['3'][index])
                cnt_fulldict['4'].append(cnt_dict['4'][index])
                cnt_fulldict['5'].append(cnt_dict['5'][index])
                cnt_fulldict['6'].append(cnt_dict['6'][index])
                break

        if not found:
            bc_fulldict['1'].append('-')
            bc_fulldict['2'].append('-')
            bc_fulldict['3'].append('-')
            bc_fulldict['4'].append('-')
            bc_fulldict['5'].append('-')
            bc_fulldict['6'].append('-')

            cnt_fulldict['1'].append(0)
            cnt_fulldict['2'].append(0)
            cnt_fulldict['3'].append(0)
            cnt_fulldict['4'].append(0)
            cnt_fulldict['5'].append(0)
            cnt_fulldict['6'].append(0)

    # adds the barcode information to the loom file
    ds.ca['linBarcode_1'] = np.array(bc_fulldict['1'], dtype='S%r' % len_bc)
    ds.ca['linBarcode_2'] = np.array(bc_fulldict['2'], dtype='S%r' % len_bc)
    ds.ca['linBarcode_3'] = np.array(bc_fulldict['3'], dtype='S%r' % len_bc)
    ds.ca['linBarcode_4'] = np.array(bc_fulldict['4'], dtype='S%r' % len_bc)
    ds.ca['linBarcode_5'] = np.array(bc_fulldict['5'], dtype='S%r' % len_bc)
    ds.ca['linBarcode_6'] = np.array(bc_fulldict['6'], dtype='S%r' % len_bc)

    # adds the count information to the loom file
    ds.ca['linBarcode_count_1'] = np.array(cnt_fulldict['1'], dtype=int)
    ds.ca['linBarcode_count_2'] = np.array(cnt_fulldict['2'], dtype=int)
    ds.ca['linBarcode_count_3'] = np.array(cnt_fulldict['3'], dtype=int)
    ds.ca['linBarcode_count_4'] = np.array(cnt_fulldict['4'], dtype=int)
    ds.ca['linBarcode_count_5'] = np.array(cnt_fulldict['5'], dtype=int)
    ds.ca['linBarcode_count_6'] = np.array(cnt_fulldict['6'], dtype=int)

    ds.close()


def log(*args, **kwargs):
    import sys
    print(*args, **kwargs, file=sys.stderr)


def main():
    args = parse_arguments()

    input_dir = args.path
    output_dir = args.name
    barcode_start = args.start
    barcode_end = args.end
    len_bc = barcode_end - barcode_start
    minham = args.hamming + 1

    os.makedirs(output_dir)

    # PART I + II: Barcode extraction and reads construction

    # 1. Extracts reads aligning to barcode-chromosome,
    # 2. extracts barcodes, UMIs and cellIDs from reads,
    # 3. outputs UMI-sorted reads with barcodes

    cell_ids = read_cellid_barcodes(
        os.path.join(input_dir, 'filtered_gene_bc_matrices', args.genome_name, 'barcodes.tsv'))

    sorted_reads = read_bam(
        os.path.join(input_dir, 'possorted_genome_bam.bam'),
        os.path.join(output_dir, args.chromosome + '_entries.bam'),
        cell_ids, args.chromosome, barcode_start, barcode_end)

    log(f'Read {len(sorted_reads)} reads containing (parts of) the barcode')
    with open(os.path.join(output_dir, 'reads.txt'), 'w') as read_file:
        print(
            '#Each output line corresponds to one read and has the following style: '
            'CellID\tUMI\tBarcode\n'
            '# dash (-) = barcode base outside of read, '
            '0 = deletion in barcode sequence (position unknown)', file=read_file)
        for read in sorted_reads:
            print(read.cell_id, read.umi, read.barcode, sep='\t', file=read_file)

    # Part III: Molecule construction

    # 1. Forms groups of reads with identical CellIDs and UMIs => belong to one molecule,
    # 2. forms consensus sequence of all barcodes of one group,
    # 3. outputs molecules and corresponding CellIDs/UMIs

    molecules = compute_molecules(sorted_reads)
    log(f'Detected {len(molecules)} molecules')
    with open(os.path.join(output_dir, 'molecules.txt'), 'w') as mol_file:
        print(
            '#Each output line corresponds to one molecule and has the following style: '
            'CellID\tUMI\tBarcode\n'
            '# dash (-) = barcode base outside of read, '
            '0 = deletion in barcode sequence (position unknown)', file=mol_file)
        for molecule in molecules:
            print(molecule.cell_id, molecule.umi, molecule.barcode, sep='\t', file=mol_file)

    # Part IV: Cell construction

    # 1. Forms groups of molecules (with set barcode minimum length) that have identical cellIDs
    #    => belong to one cell,
    # 2. counts number of appearances of each barcode in each group,
    # 3. starting from the barcode with the lowest count, compares to barcodes starting with
    #    the highest counts of a group and calculates hamming distance. If distance is below
    #    threshold, the two barcodes and counts are merged. Repetition until no barcodes with
    #    hamming distance below threshold can be found (note that this way of merging is greedy),
    # 4. Outputs for each cells all its barcodes and corresponding counts

    cells = compute_cells(molecules, args.min_length, minham)
    log(f'Detected {len(cells)} cells')
    with open(os.path.join(output_dir, 'cells.txt'), 'w') as cell_file:
        print(
            '#Each output line corresponds to one cell and has the following style: '
            'CellID\tBarcode1\tCount1\tBarcode2\tCount2...\n'
            '# dash (-) = barcode base outside of read, '
            '0 = deletion in barcode sequence (position unknown)', file=cell_file)
        for cell in cells:
            row = [cell.cell_id, ':']
            sorted_barcodes = sorted(cell.barcode_counts, key=lambda x: cell.barcode_counts[x], reverse=True)
            for barcode in sorted_barcodes:
                row.extend([barcode, cell.barcode_counts[barcode]])
            print(*row, sep='\t', file=cell_file)

    # Part V + VI: Barcodes filtering and grouping

    cells = filter_cells(cells, molecules)
    with open(os.path.join(output_dir, 'cells_filtered.txt'), 'w') as cellfilt_file:
        print(
            '#Each output line corresponds to one cell and has the following style: '
            'CellID\t:\tBarcode1\tCount1\tBarcode2\tCount2...\n'
            '# dash (-) = barcode base outside of read, '
            '0 = deletion in barcode sequence (position unknown)', file=cellfilt_file)

        for cell in cells:
            sort_d = sorted(cell.barcode_counts.items(), key=operator.itemgetter(1), reverse=True)
            if not sort_d:
                continue
            row = [cell.cell_id, ':']
            for tup in sort_d:
                row.extend(tup)
            print(*row, sep='\t', file=cellfilt_file)
    # cellIDs and filtered barcodes can be found in cells_filtered.txt

    # 2. Groups cells with same barcodes that most likely stem from one clone. Outputs a file
    #    with all clones and cellIDs belonging to each clone

    groups_dict = defaultdict(list)
    for cell in cells:
        # forms groups of cells with same barcode
        for barcode, count in cell.barcode_counts.items():
            groups_dict[barcode].append((cell.cell_id, count))

    log(f'Detected {len(groups_dict)} unique cell groups')
    # in groups.txt all barcodes and their corresponding cellIDs can be found
    with open(os.path.join(output_dir, 'groups.txt'), 'w') as groups_file:
        print(
            '#Each output line corresponds to one barcode group (clone) and has '
            'the following style: Barcode\t:\tCellID1\tbarcode-count1\tCellID2\tbarcode-count2...\n'
            '# dash (-) = barcode base outside of read, '
            '0 = deletion in barcode sequence (position unknown)', file=groups_file)

        for barcode in sorted(groups_dict):
            row = [barcode, ':']
            for cell_id, barcode_count in groups_dict[barcode]:
                row.extend([cell_id, barcode_count])
            print(*row, sep='\t', file=groups_file)

    # Create a loom file if requested
    if args.loom:
        # TODO temporary
        cell_col = []
        for cell in cells:
            cell_col.append(cell.cell_id)
            cell_col.append(cell.barcode_counts)
        write_loom(cell_col, input_dir, output_dir, len_bc)

    log('Run completed!')


if __name__ == '__main__':
    main()

"""
Description: Program for the extraction and filtering of random barcodes from single-cell
sequencing data

Preparation: Program processes cell ranger output files. Run cell ranger before. Follow
instructions => cellranger_instructions.sh

Run: Run program in cellranger 'outs' directory OR indicate path to 'outs'-directory via --path flag
"""

import pysam
import os.path
import argparse
import numpy as np
from collections import Counter
from collections import defaultdict
import operator


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


def read_bam(bam_path, output_bam_path, cell_ids, chr_name, barcode_start, barcode_end):
    """
    bam_path -- path to input BAM file
    output_bam_path -- path to an output BAM file. All reads on the chromosome that have the
        required tags are written to this file
    """
    alignment_file = pysam.AlignmentFile(bam_path)
    out_bam = pysam.AlignmentFile(output_bam_path, 'wb', template=alignment_file)

    # Fetches those reads aligning to the artifical, barcode-containing chromosome
    read_col = []
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
        umi = read.get_tag('UB')
        read_col.append((cell_id, umi, barcode))

    # sorts reads first based on UMI, then CellID, then barcode
    sorted_reads = sorted(read_col, key=lambda read: (read[1], read[0], read[2]))
    alignment_file.close()
    out_bam.close()
    assert len(sorted_reads) == 0 or len(sorted_reads[0][2]) == barcode_end - barcode_start
    return sorted_reads


def compute_molecules(read_sorted):
    if read_sorted:
        len_bc = len(read_sorted[0][2])
    else:
        len_bc = -1
    # extracts the start and end index of groups with identical UMI and cellID
    group_pos = [0]
    for i in range(0, len(read_sorted) - 1):
        if not (read_sorted[i][1] == read_sorted[i + 1][1] and read_sorted[i][0] == read_sorted[i + 1][0]):
            group_pos.append(i + 1)

    # creates a list of sublists, each representing one group of reads with identical UMI/cellID
    groups = []
    for i in range(0, len(group_pos) - 1):
        groups.append(read_sorted[group_pos[i]:group_pos[i + 1]])
    groups.append(read_sorted[group_pos[-1]:(len(read_sorted) + 1)])

    # converts each sequence of each group that is greater than 1 into a binary code, sums up binary code of all sequences, calculates max value for each position and outputs consensus sequence
    letters = np.array(['A', 'C', 'G', 'T', '-', '0'])

    mol_col = list()
    for group in groups:  # takes out each group
        if len(group) > 1:  # filters out groups that contain only one read
            consens_np = np.zeros([len_bc, 6], dtype='float16')
            for j in range(0, len(group)):  # takes out each sequence from a group
                align = np.zeros([len_bc, 6], dtype='float16')
                for (l, s) in enumerate(group[j][2]):  # takes out each base from sequence
                    # turns each base into a number and position in numpy array
                    if s == 'A':
                        align[l, 0] = 1
                    elif s == 'C':
                        align[l, 1] = 1
                    elif s == 'G':
                        align[l, 2] = 1
                    elif s == 'T':
                        align[l, 3] = 1
                    elif s == '-':
                        align[l, 4] = 0.1
                    elif s == '0':
                        align[l, 5] = 0.1
                consens_np = consens_np + align  # sums up numbers of each position
            # calculate base with maximum count for each position
            bin_consens = np.argmax(align, axis=1)
            # converts maximum counts into consensus sequence
            x = letters[bin_consens]
            consensus = ''.join(x)

            mol_col.append((group[0][0], group[0][1], consensus))
        else:
            mol_col.append((group[0][0], group[0][1], group[0][2]))
    # calling mol_col will give a list of all molecules with corresponding UMIs/cellIDs. See molecules.txt

    # sorts molecules based on cellIDs, then barcodes, then UMIs
    mol_sorted = sorted(mol_col, key=lambda mol: (mol[0], mol[2], mol[1]))
    return groups, mol_sorted


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
    import loompy

    loom_name = os.path.basename(input_dir[:-5])
    pwd_loom = os.path.join(run_name, loom_name + '.loom')
    if not os.path.exists(pwd_loom):
        loompy.create_from_cellranger(input_dir[:-5], run_name)
    # connects to the just created loom file in order to modify it
    ds = loompy.connect(pwd_loom)
    # gets a list of all cellIDs appearing in the loom file
    all_cellIDs = ds.ca.CellID

    # brings barcode data into correct format for loom file. Array must have same shape as all_cellIDs
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


def main():
    args = parse_arguments()

    chr_name = args.chromosome
    input_dir = args.path
    output_dir = args.name
    start_bc = args.start
    end_bc = args.end
    len_bc = end_bc - start_bc
    minlen_bc = args.min_length - 1
    minham = args.hamming + 1

    # Creating an output folder named after user or default defined run-name in current working directory
    os.makedirs(output_dir)

    # PART I + II: Barcode extraction and reads construction

    # 1. Extracts reads aligning to barcode-chromosome,
    # 2. extracts barcodes, UMIs and cellIDs from reads,
    # 3. outputs UMI-sorted reads with barcodes

    cell_ids = read_cellid_barcodes(os.path.join(input_dir, 'filtered_gene_bc_matrices', args.genome_name, 'barcodes.tsv'))

    sorted_reads = read_bam(
        os.path.join(input_dir, 'possorted_genome_bam.bam'),
        os.path.join(output_dir, chr_name + '_entries.bam'),
        cell_ids, chr_name, start_bc, end_bc)

    with open(os.path.join(output_dir, 'reads.txt'), 'w') as read_file:
        read_file.write(
            '#Each output line corresponds to one read and has the following style: CellID\tUMI\tBarcode' + '\n' + '# dash (-) = barcode base outside of read, 0 = deletion in barcode sequence (position unknown)' + '\n')
        for read in sorted_reads:
            print(*read[:3], sep='\t', file=read_file)

    # Part III: Molecule construction

    # 1. Forms groups of reads with identical CellIDs and UMIs => belong to one molecule,
    # 2. forms consensus sequence of all barcodes of one group,
    # 3. outputs molecules and corresponding CellIDs/UMIs

    groups, mol_sorted = compute_molecules(sorted_reads)
    with open(os.path.join(output_dir, 'molecules.txt'), 'w') as mol_file:
        mol_file.write(
            '#Each output line corresponds to one molecule and has the following style: CellID\tUMI\tBarcode' + '\n' + '# dash (-) = barcode base outside of read, 0 = deletion in barcode sequence (position unknown)' + '\n')
        for mol in mol_sorted:
            print(*mol[:3], sep='\t', file=mol_file)

    # Part IV: Cell construction

    # 1. Forms groups of molecules (with set barcode minimum length) that have identical cellIDs
    #    => belong to one cell,
    # 2. counts number of appearances of each barcode in each group,
    # 3. starting from the barcode with the lowest count, compares to barcodes starting with
    #    the highest counts of a group and calculates hamming distance. If distance is below
    #    threshold, the two barcodes and counts are merged. Repetition until no barcodes with
    #    hamming distance below threshold can be found (note that this way of merging is greedy),
    # 4. Outputs for each cells all its barcodes and corresponding counts

    barcode_list = []
    cellid_list = []

    for mol in mol_sorted:
        cellid = mol[0]
        umi = mol[1]
        barcode = mol[2]
        pure_bc = barcode.strip('-')
        pure_bc0 = barcode.strip('0')
        if len(pure_bc) > minlen_bc and len(
                pure_bc0) > minlen_bc:  # filters out barcodes shorter than min length
            barcode_list.append(barcode)
            cellid_list.append(cellid)

    # extracts the start and end index of groups with identical cellID
    group_pos = [0]
    for i in range(0, len(cellid_list) - 1):
        if cellid_list[i] != cellid_list[i + 1]:
            group_pos.append(i + 1)

    # creates a list of sublists, each representing one group of molecules with identical cellID
    cellid_grp = []
    barcode_grp = []
    for i in range(0, len(group_pos) - 1):
        cellid_grp.append(cellid_list[group_pos[i]:group_pos[i + 1]])
        barcode_grp.append(barcode_list[group_pos[i]:group_pos[i + 1]])
    cellid_grp.append(cellid_list[group_pos[-1]:(len(cellid_list) + 1)])
    barcode_grp.append(barcode_list[group_pos[-1]:(len(barcode_list) + 1)])

    cell_file = open(os.path.join(output_dir, 'cells.txt'), 'w')
    cell_file.write(
        '#Each output line corresponds to one cell and has the following style: CellID\tBarcode1\tCount1\tBarcode2\tCount2...' + '\n' + '# dash (-) = barcode base outside of read, 0 = deletion in barcode sequence (position unknown)' + '\n')

    # merges barcodes and counts below hamming distance
    cell_col = []
    found = False
    for group in cellid_grp:
        cellid = group[0]
        bcgrp = barcode_grp[cellid_grp.index(group)]
        bc_counts = Counter(bcgrp)  # counts the appearances of different barcodes in each group
        results = defaultdict(int)
        mc = sorted(bc_counts.most_common(), key=lambda x: -len(x[0].strip('-')))  # sorts barcodes based on counts
        while True:
            x, n = mc.pop(-1)  # takes out and remove barcode with lowest count from list
            if len(mc) == 0:  # or '0' in x: #if barcode is the last in the list or it contains insertions/deletions (cannot be compared) just keeps barcode without merging
                results[x] += n
                break
            for i, m in mc:  # goes through remaining barcodes in list
                hamming = 0
                overlap_count = 0
                for l, k in zip(x, i):
                    if l != '-' and k != '-':  # only compares base-containing and not empty position
                        overlap_count += 1  # counts the overlap of two barcodes
                        if l != k:
                            hamming += 1  # calculates hamming distance based on the similarity of each base-pair
                if hamming < minham and overlap_count != 0:  # filters out barcode-pairs with hamming distance below set threshold or with no overlapp
                    if len(i.strip('-')) == len_bc:  # only full barcodes are merged with other groups
                        results[i] += n
                    else:
                        results[x] += n
                    found = True
                    break

            if not found:  # barcodes that never undergo the hamming distance threshold, are not merged
                results[x] += n
            else:
                found = False

        cell_col.append(cellid)
        cell_col.append(results)
        cell_file.write(cellid + '\t:\t')
        results_sorted = sorted(results, key=lambda x: -results[x])
        for key in results_sorted:
            cell_file.write(key + '\t' + str(results[key]) + '\t')
        cell_file.write('\n')
    cell_file.close()

    # calling cell_col will give all cells and filtered barcodes. See below.


    # Part V + VI: Barcodes filtering and grouping

    # 1. Filters barcodes according to two criteria:
    #    a) Barcodes that have only a count of one and can be found in another cell are most
    #       likely results of contamination and are removed,
    #    b) Barcodes that have only a count of one and are also only based on one read are
    #       also removed
    # 2. Groups cells with same barcodes that most likely stem from one clone. Outputs a file
    #    with all clones and cellIDs belonging to each clone

    bc_all = list()
    cellids = list()
    for i in range(1, len(cell_col), 2):
        for key in cell_col[i]:
            cellids.append(cell_col[i - 1])
            bc_all.append(key)
            bc_all.append(cell_col[i][key])

    # filters out barcodes with a count of one that appear in another cell
    for i in range(1, len(cell_col), 2):
        dict_cp = dict(cell_col[i])
        for key in dict_cp:
            if key in cell_col[i]:
                if cell_col[i][key] == 1:  # filters out barcodes that are only based on one molecule
                    if bc_all.count(
                            key) > 1:  # filters out barcodes that appear more than once in the whole list
                        del cell_col[i][key]  # removes barcodes that meet both criteria
                    else:
                        for j in groups:  # groups is a list of groups of reads with identical UMIs/cellIDs (see part II)
                            if (cell_col[i - 1] == j[0][
                                0]):  # if cellID is identical to cellID in groups, it keeps the group
                                if len(
                                        j) == 1:  # filters out those barcodes that are based on only one read => group has only a length of one
                                    if key in cell_col[i]:
                                        del cell_col[i][key]  # deletes those barcodes

    # calling cell_col will give a list of all cellIDs and only the filtered barcodes

    with open(os.path.join(output_dir, 'cells_filtered.txt'), 'w') as cellfilt_file:
        cellfilt_file.write(
            '#Each output line corresponds to one cell and has the following style: CellID\t:\tBarcode1\tCount1\tBarcode2\tCount2...' + '\n' + '# dash (-) = barcode base outside of read, 0 = deletion in barcode sequence (position unknown)' + '\n')

        groups_dict = dict()
        for i in range(0, len(cell_col), 2):
            sort_d = sorted(cell_col[i + 1].items(), key=operator.itemgetter(1))
            sort_d.reverse()
            if len(sort_d) != 0:
                cellfilt_file.write(cell_col[i] + '\t:\t')
                for tup in sort_d:
                    cellfilt_file.write(tup[0] + '\t' + str(tup[1]) + '\t')
                cellfilt_file.write('\n')
            # cellIDs and filtered barcodes can be found in cells_filtered.txt

            # forms groups of cells with same barcode
            for key in cell_col[i + 1]:
                if not key in groups_dict:  # creates a new group if not existing yet. Saves cellID in a list
                    groups_dict.update({key: [cell_col[i], cell_col[i + 1][key]]})
                else:  # updates an existing group by appending cellID to the cellID list
                    groups_dict[key].append(cell_col[i])
                    groups_dict[key].append(cell_col[i + 1][key])

    groupsdict_s = sorted(groups_dict.items(), key=operator.itemgetter(0))

    # in groups.txt all barcodes and their corresponding cellIDs can be found
    with open(os.path.join(output_dir, 'groups.txt'), 'w') as groups_file:
        groups_file.write(
            '#Each output line corresponds to one barcode group (clone) and has the following style: Barcode\t:\tCellID1\tbarcode-count1\tCellID2\tbarcode-count2...' + '\n' + '# dash (-) = barcode base outside of read, 0 = deletion in barcode sequence (position unknown)' + '\n')

        for i in range(0, len(groupsdict_s)):
            groups_file.write(groupsdict_s[i][0] + '\t:\t')
            for j in range(0, len(groupsdict_s[i][1]), 2):
                groups_file.write(groupsdict_s[i][1][j] + '\t' + str(groupsdict_s[i][1][j + 1]) + '\t')
            groups_file.write('\n')

    # Create a loom file if requested
    if args.loom:
        write_loom(cell_col, input_dir, output_dir, len_bc)

    print('Run completed!')


if __name__ == '__main__':
    main()

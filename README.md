![CI status](https://github.com/frisen-lab/TREX/workflows/CI/badge.svg)

# Overview

TREX is an experimental workflow that enables simultaneous lineage TRacking and EXpression profiling of single cells using RNA-sequencing. 

An essential part of this workflow is presented here: the extraction of genetic barcodes or "cloneIDs" from single-cell transcriptomes and the reconstruction of related cells.

The computational tool uses BAM files of one or multiple sequencing libraries as an input for the generation of cloneID count matrices and identifies clonally related cells based on Jaccard similarity between each pair of cloneID+cells. 

Currently, TREX is compatible with common RNA-sequencing library preparation methods and data formats provided by 10X Chromium and 10X Visium.


# Installation

TREX requires Python 3.7 or newer.
We recommend that you install TREX into a "virtual environment", which can be done by running the following commands:

    python3 -m venv trex-venv
    trex-venv/bin/pip install git+https://github.com/frisen-lab/TREX.git

`trex-venv` is the name of the directory that will be created and which will contain the virtual environment.
You can choose a different name.

Activate the virtual environment by running
`source trex-venv/bin/activate`.
You need to repeat this step in every new shell in order to be able to run TREX.
Finally, test the installation by running `trex --version`.


# Running TREX on a minimal test dataset

Clone the Git repository or [download it as a ZIP file](https://github.com/frisen-lab/TREX/archive/main.zip) and unpack it.
The directory `tests/data/` contains a test dataset in the proper format.
Run TREX on it:

    trex run10x -s 695 -e 724 tests/data/

This will create a folder `trex_run/` in the current directory
(use the `-o` option to choose a different folder) with the results.

See the "Runnning TREX" section below for further details.


# Reproducing results from the Nature Neuroscience paper

See the manuscript at <https://doi.org/10.1038/s41593-022-01011-x>.
This section lists some details not available in the paper.
See the `references/` directory for the necessary files.

The FASTA reference used with CellRanger was created by appending
`H2b-EGFP-30N-LTR.fa` to the end of the GRCm38 (mm10) reference FASTA:

    cat genome.fa chrH2B-EGFP-N.fa > mm10_H2B-EGFP-30N_genome.fa

The GTF annotations file was created by appending `chrH2B-EGFP-N.gtf`
to the existing annotations:

    cat genes.gtf chrH2B-EGFP-N.gtf > mm10_H2B-EGFP-30N_genes.gtf

Then a new CellRanger reference can be created:

    cellranger mkref --genome=mm10_H2B-EGFP-30N --fasta=mm10_H2B-EGFP-30N_genome.fa --genes=mm10_H2B-EGFP-30N_genes.gtf > mkref_mm10_H2B-EGFP-30N.out


# Running TREX

The input directory for TREX must be a Cell Ranger output directory.
See the contents of the `tests/data/outs` directory to learn which are the minimum files necessary.
Cell Ranger must have been configured to map against a reference augmented by an extra chromosome that contains the cloneID. By default, that extra chromosome is assumed to be the last in the BAM file (use `--chromosome` to choose a different one).
The options `-s` and `-e` set where on the extra chromosome the cloneID is located (`-s` gives start and `-e` gives end in 1-based coordinates).

Please also run `trex run10x --help` to see the other available command-line options.


## Pipeline steps overview

This is an overview of the steps that the `trex run10x` command performs.

1. Retrieve usable reads from the input BAM file.
   A usable read fulfills these requirements:
   - It aligns to the region specified with the `--chr`, `-s` and
     `-e` flags or, if the flags are not given, to the region that
     has been automatically identified to be the region containing
     the variable CloneID sequence,
   - it has both an associated cell ID and UMI (SAM tags `CB` and `UB`),
   - its cell ID is included in the list of allowed cell IDs
     (if such a list is provided with `--filter-cellid` or `-f`).
2. Group reads with identical cell ID and UMI into *molecules*.
   The clone ID of the molecule is taken to be the consensus of the clone IDs.
3. Error-correct clone IDs of the molecules.
4. Group molecules by cell ID into cells.
5. Filter bad clone IDs in each cell:
   Rare clone IDs also found in another cell are considered to be contaminants.
   Clone IDs supported by only a single read are removed.
6. Cluster the cells into clones by creating a *clone graph*.
   Edges are drawn between cells that appear to belong to the same clone.
   The connected components of the graph are considered to be the clones.
   (A clone is thus simply a set of cells.)
7. Error-correct the clone graph by removing spurious edges ("bridges").


## Output files

`trex run10x` by default writes its results to a newly created directory
named `trex_run`.
The name of the output directory can be changed with `--output` (or `-o`).
The files created in the output directory are described below.
They are listed in the order in which the program creates them,
see the "Pipeline steps overview" section for more details about what is
done in each step.

Many of the files are tables in either tab-separated value (TSV) format
or in comma-separated values (CSV) format.


### `log.txt`

This file contains a copy of the output that `trex run10x` prints to the
terminal.

### `entries.bam`

All usable reads from the input BAM file.
(See the pipeline steps overview for a description of "usable reads".)


### `reads.txt`

A table listing *cell_id*, *umi* and *clone_id* for each usable input read.
Example:

    #cell_id          umi         clone_id
    TGACGGCGTTACCAGT  AAAAAACTGT  TGTCAATCGTTCGGTTGAGCAAGATCTTAG


### `molecules.txt`

A table listing *cell_id*, *umi* and *clone_id* for each detected molecule.
Example:

    #cell_id                umi             clone_id
    AAACCTGAGAGGTACC        AGTTAAAGTA      TGTCAATCGTTCGGTTGAGCAAGATCTTAG


### `molecules_corrected.txt`

A table with molecules where clone IDs have been error corrected.

The file is the same as `molecules.txt` with some changes in the *clone_id*
column.


### `cells.txt`

A table listing the detected cells.
The fields are *cell_id*, a colon (`:`) and then
pairs of columns *clone_id1* and *count1* for each clone ID found in that cell.
Example:

    #cell_id          :  clone_id1                       count1  clone_id2 count2  ...
    AAACCTGAGAGGTACC  :  TGTCAATCGTTCGGTTGAGCAAGATCTTAG  1
    AAACCTGAGTAGCGGT  :  TGTCAATCGTTCGGTTGAGCAAGATCTTAG  3


### `cells_filtered.txt`

The same as `cells.txt`, but with error-corrected clone ID counts.
Cells that end up without any clone IDs after error correction are removed.


### `umi_count_matrix.csv`

A matrix of UMI counts with cells as columns and clone IDs as rows.

Only created if option `--umi-matrix` is used.


### `components.txt`

Connected components of the clone graph.


### `components_corrected.txt`

Connected components of the clone graph after error correction.


### `graph.pdf`, `graph_corrected.pdf`

Plots of the clone graph and the error-corrected clone graph.

These are only created if option `--plot` is used.

`graph_corrected.gv` and `graph.gv` are textual descriptions of the graphs in GraphViz format.


### `clones.txt`

A table listing the cell IDs belonging to each clone.
The columns are *clone_id* and *cell_id* where *clone_id* is a number that
identifies the clone.

    #clone_id,cell_id
    1,TGGCGCAAGAATAGGG


### `clone_sequences.txt`

A table listing the 30N sequence of each clone ID.
The columns are *clone_id* and *clone_seq* where *clone_id* is a number
that identifies the clone and *clone_seq* its nucleotide sequence.

    #clone_id,clone_seq
    1,ACTAGGAGATTGACGGATCACCTTTGGTCG


### `data.loom`

A [loom file](http://linnarssonlab.org/loompy/).

This file is created only if option `--loom` (or `-l`) is used.

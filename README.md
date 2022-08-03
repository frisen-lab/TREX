![CI status](https://github.com/frisen-lab/TREX/workflows/CI/badge.svg)

# Overview

TREX is an experimental workflow that enables simultaneous lineage TRacking and EXpression profiling of single cells using RNA-sequencing. The method is described in the paper [Clonal relations in the mouse brain revealed by single-cell and spatial transcriptomics](https://doi.org/10.1038/s41593-022-01011-x).

An essential part of this workflow is presented here: the extraction of genetic barcodes or "cloneIDs" from single-cell transcriptomes and the reconstruction of related cells.

The tool uses BAM files of one or multiple sequencing libraries as an input for the generation of cloneID count matrices and identifies clonally related cells based on Jaccard similarity between each pair of cloneID+cells.

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


# Running TREX on data from the Nature Neuroscience paper

Here we show how to run TREX on the data we analyzed in our Nature Neuroscience paper (<https://doi.org/10.1038/s41593-022-01011-x>).

These instructions have been tested with Cell Ranger 6.1.2, but the original analysis was done with Cell Ranger 2.2.0.


## Data

The data is available under GEO accession [GSE153424](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153424). The instructions below analyze sample [GSM4644060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4644060) because it is relatively small.

GSM4644060 is also called "brain1_str" (*str* stands for striatum) in the GEO description, and it has the ID "10X_41" in [Supplementary Table 4](https://static-content.springer.com/esm/art%3A10.1038%2Fs41593-022-01011-x/MediaObjects/41593_2022_1011_MOESM6_ESM.xlsx) in the paper.

As can be seen on the [overview page for GSM4644060](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4644060), reads for the dataset are available from [SRA experiment accession SRX8627776](https://www.ncbi.nlm.nih.gov/sra?term=SRX8627776), which in turn links to two run accessions:
* [SRR12103475](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR12103475)
* [SRR12103476](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR12103476)

We use the run accessions to retrieve the data.


## Prerequisites

- Go to the [Cell Ranger download page](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest). Download and install Cell Ranger.
- Download and extract the mouse reference dataset `refdata-gex-mm10-2020-A.tar.gz`:

       curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
       tar xvf refdata-gex-mm10-2020-A.tar.gz

- Create a custom reference by adding an extra chrH2B-EGFP-N contig (EGFP-cloneID virus) to the above mouse reference dataset.
  `references/` refers to the directory in this repository. `cellranger mkref` takes about one hour.
  You may choose to continue to the next step (downloading reads) while this is running.

      cat refdata-gex-mm10-2020-A/fasta/genome.fa references/chrH2B-EGFP-N.fa > mm10_H2B-EGFP-30N_genome.fa
      cat refdata-gex-mm10-2020-A/genes/genes.gtf references/chrH2B-EGFP-N.gtf > mm10_H2B-EGFP-30N_genes.gtf
      cellranger mkref --genome=mm10_H2B-EGFP-30N --fasta=mm10_H2B-EGFP-30N_genome.fa --genes=mm10_H2B-EGFP-30N_genes.gtf > mkref_mm10_H2B-EGFP-30N.out


## Downloads reads

Create a `fastq` directory and change into it:

    mkdir fastq
    cd fastq

If you do not have it, install `fastq-dump`. For example, if you have [Conda](https://docs.conda.io/) (with the [Bioconda channels activated](http://bioconda.github.io/user/install.html#set-up-channels)), run

    conda create -n trex sra-tools
    conda activate trex

Download the reads:

    fastq-dump --gzip --split-3 --defline-qual '+' --defline-seq '@$ac.$sn' SRR12103475 SRR12103476

Rename the files so that Cell Ranger can find them:

    mv SRR12103475_1.fastq.gz SRR12103475_S1_L001_R1_001.fastq.gz
    mv SRR12103475_2.fastq.gz SRR12103475_S1_L001_R2_001.fastq.gz
    mv SRR12103476_1.fastq.gz SRR12103476_S1_L001_R1_001.fastq.gz
    mv SRR12103476_2.fastq.gz SRR12103476_S1_L001_R2_001.fastq.gz

    cd ..


## Run Cell Ranger

Continue with this step only when the above steps have finished.
Run `cellranger count`:

    cellranger count --transcriptome=mm10_H2B-EGFP-30N --id=brain1_str --fastqs=fastq/ --sample=SRR12103475 --sample=SRR12103476 --expect-cells=2299

The `--expect-cells` parameter is set to the number of cells loaded per well in the 10X chip used for droplet generation (4000 in this case, see the "cells recovered for 10X" column in Suppl. Table 4) divided by 1.74.


## Run TREX

With `trex` installed (as described above), run

    trex run10x -o trex_brain1_str brain1_str

Results will be written to a new directory named `trex_brain1_str`.


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
     the variable cloneID sequence,
   - it has both an associated cell ID and UMI (SAM tags `CB` and `UB`),
   - its cell ID is included in the list of allowed cell IDs
     (if such a list is provided with `--filter-cellid` or `-f`).
2. Group reads with identical cell ID and UMI into *molecules*.
   The cloneID of the molecule is taken to be the consensus of the cloneIDs.
3. Error-correct cloneIDs of the molecules.
4. Group molecules by cell ID into cells.
5. Filter bad cloneIDs in each cell:
   Rare cloneIDs also found in another cell are considered to be contaminants.
   CloneIDs supported by only a single read are removed.
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

A table with molecules where cloneIDs have been error corrected.

The file is the same as `molecules.txt` with some changes in the *clone_id*
column.


### `cells.txt`

A table listing the detected cells.
The fields are *cell_id*, a colon (`:`) and then
pairs of columns *clone_id1* and *count1* for each cloneID found in that cell.
Example:

    #cell_id          :  clone_id1                       count1  clone_id2 count2  ...
    AAACCTGAGAGGTACC  :  TGTCAATCGTTCGGTTGAGCAAGATCTTAG  1
    AAACCTGAGTAGCGGT  :  TGTCAATCGTTCGGTTGAGCAAGATCTTAG  3


### `cells_filtered.txt`

The same as `cells.txt`, but with error-corrected cloneID counts.
Cells that end up without any cloneIDs after error correction are removed.


### `umi_count_matrix.csv`

A matrix of UMI counts with cells as columns and cloneIDs as rows.
This is a different representation of the data written to `cells_filtered.txt`.

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

A table listing the 30N sequence of each cloneID.
The columns are *clone_id* and *clone_seq* where *clone_id* is a number
that identifies the clone and *clone_seq* its nucleotide sequence.

    #clone_id,clone_seq
    1,ACTAGGAGATTGACGGATCACCTTTGGTCG


### `data.loom`

A [loom file](http://linnarssonlab.org/loompy/).

This file is created only if option `--loom` (or `-l`) is used.

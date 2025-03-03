![CI status](https://github.com/frisen-lab/TREX/workflows/CI/badge.svg)

# Overview

TREX is an experimental workflow that enables simultaneous lineage TRacking and EXpression profiling of single cells using RNA-sequencing. The method is described in the paper [Clonal relations in the mouse brain revealed by single-cell and spatial transcriptomics](https://doi.org/10.1038/s41593-022-01011-x).

An essential part of this workflow is presented here: the extraction of genetic barcodes or "cloneIDs" from single-cell or spatial transcriptomes and the reconstruction of related cells/spots.

The tool uses BAM files of one or multiple sequencing libraries as an input for the generation of cloneID count matrices and identifies clonally related cells based on Jaccard similarity between each pair of cloneID+cells.

Currently, TREX is compatible with common RNA-sequencing library preparation methods and data formats provided by 10X Chromium, 10X Visium, Smart-seq2 and 3.


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


# Changelog

See [Changelog](CHANGES.md).


# Running TREX on a minimal test dataset

Clone the Git repository or [download it as a ZIP file](https://github.com/frisen-lab/TREX/archive/main.zip) and unpack it.
The directory `tests/data/` contains a test dataset in the proper format.
Run TREX on it:

    trex run10x -s 695 -e 724 tests/data/

This will create a folder `trex_run/` in the current directory
(use the `-o` option to choose a different folder) with the results.

See the "Runnning TREX" section below for further details.


# Running TREX on single-cell data from the Nature Neuroscience paper

Here we show how to run TREX on some of the single-cell data we analyzed in our Nature Neuroscience paper (<https://doi.org/10.1038/s41593-022-01011-x>).

These instructions have been tested with Cell Ranger 6.1.2, but the original analysis was done with Cell Ranger 2.2.0.

These instructions apply only partially to the spatial data (samples GSM4644079-GSM4644094),
see the [Visium](#visium) section below for the differences.


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


## Visium

If you want to re-run our experiments on the Visium samples (GSM4644079-GSM4644094),
you need to use Space Ranger and a different virus sequence.
We used Space Ranger 1.0.0.

You can find the virus sequence in `references/pMR526_LV-EF1a-H2B-EGFP-BC.fa` in this repository.

This is the association between the GEO sample names and our custom library
names listed in Supplementary Table 4:

GEO sample id | Name in Suppl. Tab. 4
-|-
GSM4644079 | V9
GSM4644080 | V10
GSM4644081 | V11
GSM4644082 | V12
GSM4644083 | V13
GSM4644084 | V14
GSM4644085 | V15
GSM4644086 | V16
GSM4644087 | amp56
GSM4644088 | amp57
GSM4644089 | amp58
GSM4644090 | amp59
GSM4644091 | amp60
GSM4644092 | amp61
GSM4644093 | amp62
GSM4644094 | amp63

V... are the gene expression libraries and amp... are the libraries targeting the cloneID.
There will be only few reads with cloneIDs in the V... libraries.


When analyzing the Space Ranger output with TREX,
you need to add `--visium` to the command.
You can also analyze the transcriptome and amplicon data at the same time.
For example:

    trex run --visium --output trex-V9 --start 1147 --end 1176 --umi-matrix --prefix --visium --samples V9 visium-data/V9 --amplicon visium-data/amp56


# Running TREX

The input directory for TREX must be a Cell Ranger output directory. In case of Smart-Seq2 / 3 data, one BAM file with all cells or a folder with one BAM file per cell is expected (see zUMIs output)
See the contents of the `tests/data/` directory to learn which are the minimum files necessary.
Cell Ranger/zUMIs must have been configured to map against a reference augmented by an extra chromosome that contains the cloneID. By default, that extra chromosome is assumed to be the last in the BAM file (use `--chromosome` to choose a different one).
The options `-s` and `-e` set where on the extra chromosome the cloneID is located (`-s` gives start and `-e` gives end in 1-based coordinates).

Please also run `trex run10x --help` (or `trex smartseq2 --help` and `trex smartseq3 --help` respectively ) to see the other available command-line options.


## Pipeline steps overview

This is an overview of the steps that the `trex run10x`/ `trex smartseq3` command performs.

1. Retrieve usable reads from the input BAM file.
   A usable read fulfills these requirements:
   - It aligns to the region specified with the `--chr`, `-s` and
     `-e` flags or, if the flags are not given, to the region that
     has been automatically identified to be the region containing
     the variable cloneID sequence,
   - it has both an associated cell ID and UMI (SAM tags `CB` and `UB`
     in case of 10x data , SAM tags `BC` and `UB` in case of Smart-seq3 data),
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
   Edges are drawn between cells that appear to belong to the same clone
   based on the set Jaccard index threshold.
   The connected components of the graph are considered to be the clones.
   (A clone is thus simply a set of cells.)
7. Error-correct the clone graph by removing spurious edges ("bridges").

`trex smartseq2` follows a similar pipeline with the following differences:
- A useable read does not require a UMI
- Reads do not get grouped into molecules and their cloneIDs are not collapsed
  and error-corrected into consensus sequences


## Input parameters

### Jaccard index threshold

The Jaccard index measures the similarity of two sample sets, in this case
the similarity of two sets of cloneIDs. It is calculated by dividing
the number of overlapping, unique cloneIDs between cell A and B by the total
number of unique cloneIDs in cell A and B. An index of 0.0 indicates no
overlapping cloneIDs and an index of 1.0 a perfect match. The Jaccard
threshold is the Jaccard index above which two cells are merged into one
clone. It can be set with the `--jaccard-threshold` flag and is 0.7
by default, meaning cell A and B are merged into one clone if they have more
than 70% of cloneIDs in common.


### Filter cellids

Tab-separated file of cell IDs to keep in the TREX run. Adding this file via
the `--filter-cellid` or `-f` option allows to focus the analysis on specific cells
and to filter out low quality cells or doublets.
Example:

```
0	CACTCGTGGTACACACTCCG
1	CACTCGTGGTACCACAAGCA
```

## Filtering cloneIDs

Text file with cloneIDs to ignore. The format is one cloneID per line.
Adding this file via the `--filter-cloneids` option allows to ignore cloneIDs that have been identified as overrepresented or artefactual during library characterization.
Example file:

```
GGTCTCCCTATACCAACAGTATCGTCTCAA
GGGTTCTGGGATATTACGTTGACTTGAGAG
```

## Output files

TREX by default writes its results to a newly created directory
named `trex_run`.
The name of the output directory can be changed with `--output` (or `-o`).
The files created in the output directory are described below.
They are listed in the order in which the program creates them,
see the "Pipeline steps overview" section for more details about what is
done in each step.

Many of the files are tables in either tab-separated value (TSV) format
or in comma-separated values (CSV) format.


### `log.txt`

This file contains a copy of the output that a TREX run prints to the
terminal.


### `entries.bam`

All usable reads from the input BAM file.
(See the pipeline steps overview for a description of "usable reads".)


### `reads.txt`

A table listing *cell_id*, *umi* and *clone_id* for each usable input read.
Example:

    #cell_id          umi         clone_id
    TGACGGCGTTACCAGT  AAAAAACTGT  TGTCAATCGTTCGGTTGAGCAAGATCTTAG

- The character `0` in a cloneIDs signals a deleted base (CIGAR operation
  `D` in the input BAM file).
- If the read does not fully cover the cloneID region, the cloneID contains
  the character `-` for each missing base.


### `molecules.txt`

A table listing *cell_id*, *umi* and *clone_id* for each detected molecule.
Example:

    #cell_id                umi             clone_id
    AAACCTGAGAGGTACC        AGTTAAAGTA      TGTCAATCGTTCGGTTGAGCAAGATCTTAG


### `molecules_corrected.txt`

Same as `molecules.txt`, but after filtering and error correction:
- CloneIDs have been error corrected.
- Molecules are removed that did not pass filtering criteria
  (such as low-complexity cloneID filtering).

For the `run10x` subcommand, this table contains an additional `original_clone_id`
column that shows how the original cloneID before correction.


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


### `read_count_matrix.csv`

Instead of UMI count matrix, running `trex smartseq2` produces matrix of read counts with cells as columns and cloneIDs as rows.
This is a different representation of the data written to `cells_filtered.txt`.

Only created if option `--read-matrix` is used.


### `components.txt`

Connected components of the clone graph.


### `components_corrected.txt`

Connected components of the clone graph after error correction.


### `graph.pdf`, `graph_corrected.pdf`

Plots of the clone graph and the error-corrected clone graph.

These are only created if option `--plot` is used.

`graph_corrected.gv` and `graph.gv` are textual descriptions of the graphs in GraphViz format.

### `doublets.txt`

A list of the cell IDs of the cells that were detected to be doublets.


### `clones.txt`

A table listing the cell IDs belonging to each clone.
The columns are *clone_id* and *cell_id* where *clone_id* is a number that
identifies the clone.

    clone#   cell_id
    1        TGGCGCAAGAATAGGG


### `clone_details.txt`

A table listing all clones (a clone is a set of cells), one line per clone.
The columns are:

* *clone_nr*: A number identifying the clone
* *clone_id*: The most common cloneID among all cells in the clone.
* *n_cells*: The number of cells in the clone.
* *clone_ids_per_cell*: The average number of cloneIDs per cell.

The columns are *clone#* and *clone_seq* where *clone_id* is a number
that identifies the clone and *clone_seq* its nucleotide sequence.

    clone_nr  clone_id                        n_cells clone_ids_per_cell
    1         ACTAGGAGATTGACGGATCACCTTTGGTCG  3       1.00


### `data.loom`

A [loom file](http://linnarssonlab.org/loompy/).

This file is created only if option `--loom` (or `-l`) is used.


# Creating a quality control report

```shell
trex qc --plot-jaccard-matrix --plot-hamming-distance DIRECTORY
```

*qc* takes as an input the directory (or directories) of trex output.
Plotting the jaccard similarity matrix between cells requires some time as jaccard similarity is calculated pairwise amongst all cells.
This can be activated adding the optional flag `--plot-jaccard-matrix`.
Hamming distance between all viral cloneIDs found in the dataset after each step can be plotted by means of the optional flag `--plot-hamming-distance`.

This will add a PDF file named *quality_report.pdf* describing the quality of the TREX run inside the same folder with the TREX output.

This report contains:

### Overall results

- Histogram of clone sizes
- Histogram of how many unique cloneIDs can pe found in each clone
- Histogram of how many unique cloneIDs can be found in each cell
- *(Optional)* A histogram of the Jaccard similarity values between cells and a matrix of the Jaccard similarity between all cells.
- Histogram of how many reads each detected viral cloneID molecule has.


### Per step results

Each of these plots has four subplots corresponding to different steps of the TREX pipeline.

- Histograms of how many nucleotides have been read in each molecule
- *(Optional)* Histograms of the Hamming distance between all the viral cloneIDs found
- Histograms of how many viral cloneID molecules have been found in each cell
- Histograms of how many molecules of each unique cloneID have been found in the dataset
- Histograms of How many unique cloneIDs per cell have been found


# TREX development

It is highly recommended that you develop TREX within a separate virtual environment:

    python3 -m venv --prompt trex .venv
    source .venv/bin/activate

Install TREX in "editable" mode:

    pip install -e .

Install `pre-commit` and install the pre-commit hooks
(these run at `git commit` time and do some checks on
the to-be-committed files):

    pip install pre-commit
    pre-commit install

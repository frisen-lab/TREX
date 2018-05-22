Installation
------------

This makes the dependencies available:

    python3 -m venv venv
    venv/bin/pip install -e .
    source venv/bin/activate


Minimal test dataset
--------------------

mini/ contains a small test dataset. The BAM file was created by downsampling
`/proj/uppstore2018019/P10306/pMR307_clone3_P10306_1001/10x_10_pMR307_clone3/outs/possorted_genome_bam.bam`:

    samtools view -s 0.01 -b 10x_10_pMR307_clone3/outs/possorted_genome_bam.bam chrEGFP-30N > new.bam

The `barcodes.tsv` is a copy of the original.

`chrEGFP-30N.fa` is a copy of `/proj/uppstore2018019/hg38_EGFP-30N/chrEGFP-30N.fa`.

Running the tool
----------------

On the minimal test dataset:

    rm -r lineage_run/
    python3 180425_lineage_loom.py -p mini -gn hg38_EGFP-30N -chr chrEGFP-30N

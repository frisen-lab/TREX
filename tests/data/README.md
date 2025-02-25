# Test dataset for run10x

The `outs/possorted_genome_bam.bam` file contains a subset of reads from an old
dataset called "pMR439_clone16". There are some issues with this data that we
describe below. Because the file still fulfills its purpose of serving as a
minimal test dataset, we have not replaced it.

The other BAM files in this subdirectory do not have the problems described
here.

## Outdated chrTomato reference

The "chrTomato" transgene sequence used in the "pMR439_clone16" experiment is
not used in current experiments. The old sequence is from a plasmid that has
been integrated into the genome using transposon-mediated integration.
In more recent experiments, a lentivirus is used.

The current chrTomato sequence can be found within this repository in the
`references/` subdirectory. To avoid confusion, the references folder no longer
contains the old, no longer relevant version of the sequence. We show it here
for reference:

```
>chrTomato-N-pMR439-outdated
GTCATCAAAGAGTTCATGCGCTTCAAGGTGCGCATGGAGGGCTCCATGAACGGCCACGAGTTCGAGATCGAGGGCGAGGG
CGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGCGGCCCCCTGCCCTTCGCCTGGGACA
TCCTGTCCCCCCAGTTCATGTACGGCTCCAAGGCGTACGTGAAGCACCCCGCCGACATCCCCGATTACAAGAAGCTGTCC
TTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGTCTGGTGACCGTGACCCAGGACTCCTCCCT
GCAGGACGGCACGCTGATCTACAAGGTGAAGATGCGCGGCACCAACTTCCCCCCCGACGGCCCCGTAATGCAGAAGAAGA
CCATGGGCTGGGAGGCCTCCACCGAGCGCCTGTACCCCCGCGACGGCGTGCTGAAGGGCGAGATCCACCAGGCCCTGAAG
CTGAAGGACGGCGGCCACTACCTGGTGGAGTTCAAGACCATCTACATGGCCAAGAAGCCCGTGCAACTGCCCGGCTACTA
CTACGTGGACACCAAGCTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAGCGCTCCGAGGGCC
GCCACCACCTGTTCCTGTACGGCATGGACGAGCTGTACAAGTAgGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGATC
TTTTTCCCTCTGCCAAAAATTATGGGGACATCATGAAGCCCCTTGAGCATCTGACTTCTGGCA
```

A GTF line would also be needed for `cellranger mkref`, which could look like
this:
```
chrTomato-N-pMR439-outdated  other  exon  1  782  .  +  .  gene_id "Tomato-N"; transcript_id "Tomato-N";
```

## Deletions every 80 nucleotides

The alignments in the BAM file contain spurious deletions every 80 nucleotides.
These exist because the pMR439_clone16 dataset was aligned against a FASTA file
where the chrTomato sequence had Windows line breaks, which the read mapper
(STAR) interpreted as an extra nucleotide. We actually considered it
advantageous to have a couple of deletions in the BAM file to test whether the
cloneID detection script can deal with them, so we left them in. (By chance,
the cloneID itself is not affected.)

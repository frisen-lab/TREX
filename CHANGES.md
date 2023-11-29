# TREX Changelog

## development version

* [#7](https://github.com/frisen-lab/TREX/issues/7):
  Set default Jaccard threshold to 0.7 so that it matches R output.
* [#33](https://github.com/frisen-lab/TREX/issues/33):
  Added a `trex qc` subcommand for quality control plots.
* [#44](https://github.com/frisen-lab/TREX/issues/44):
  Fixed some inconsistencies in the way cloneIDs are error-corrected.
* [#30](https://github.com/frisen-lab/TREX/issues/30):
  Added doublet filtering. Cells that appear to be connected to two
  subclusters are detected and removed. The cell IDs of detected doublets
  are written to `doublets.txt`. Doublet detection can be disabled with
  `--keep-doublets`.
* [#42](https://github.com/frisen-lab/TREX/issues/42):
  Added a filter to `run10x` that always removes cloneIDs with a highly
  uneven distribution of nucleotide frequencies (low-complexity cloneIDs,
  measured using Shannon entropy).
  This includes, for example, cloneIDs consisting of a single, repeated
  nucleotide such as `AAAA...`, but also `AAAAAAAAAAAAAAAAAAAAAAAAAAGAAA`.
* [#40](https://github.com/frisen-lab/TREX/issues/40):
  Added `--filter-cloneids` option for excluding certain cloneIDs
* [#26] (https://github.com/frisen-lab/TREX/issues/26):
  Changed reading of BAM files such as that it is possible to indicate a
  path with several BAM files as imput. `read_multiple` loops through all
  BAM files in a given path.
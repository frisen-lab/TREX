# TREX Changelog

## development version

* [#75](https://github.com/frisen-lab/TREX/pull/75):
  Output the GraphViz file component by component
* [#77](https://github.com/frisen-lab/TREX/pull/77):
  Do not plot clones with just two nodes
* [#79](https://github.com/frisen-lab/TREX/pull/79):
  Fix: CloneIDs where the unknown nucleotides are in the middle were not filtered
  out correctly.
* [#76](https://github.com/frisen-lab/TREX/pull/76):
  Added a `read_count` column to the `molecules*.txt` files
* `clone_sequences.txt` was renamed to `clone_details.txt`.
* Added `n_cells` an `clone_ids_per_cell` columns to `clone_details.txt`

## v0.5 (2024-08-24)

* [#30](https://github.com/frisen-lab/TREX/issues/30):
  Added doublet filtering. Cells that appear to be connected to two
  subclusters are detected and removed. The cell IDs of detected doublets
  are written to `doublets.txt`. Doublet detection can be disabled with
  `--keep-doublets`.
* [#44](https://github.com/frisen-lab/TREX/issues/44):
  Fixed some inconsistencies in the way cloneIDs are error-corrected.
* [#40](https://github.com/frisen-lab/TREX/issues/40):
  Added `--filter-cloneids` option for excluding certain cloneIDs
* [#26](https://github.com/frisen-lab/TREX/issues/26):
  Changed reading of BAM files such as that it is possible to indicate a
  path with several BAM files as imput.
* [#71](https://github.com/frisen-lab/TREX/issues/71):
  Process exclusion list (`--filter-cloneids`) a lot faster. Processing a
  long list with hundreds of thousands of entries now takes minutes instead
  of days.
* [#42](https://github.com/frisen-lab/TREX/issues/42):
  Added a filter to `run10x` that always removes cloneIDs with a highly
  uneven distribution of nucleotide frequencies (low-complexity cloneIDs,
  measured using Shannon entropy).
  This includes, for example, cloneIDs consisting of a single, repeated
  nucleotide such as `AAAA...`, but also `AAAAAAAAAAAAAAAAAAAAAAAAAAGAAA`.

## v0.4 (2023-10-10)

* [#7](https://github.com/frisen-lab/TREX/issues/7):
  Set default Jaccard threshold to 0.7 so that it matches R output.
* [#33](https://github.com/frisen-lab/TREX/issues/33):
  Added a `trex qc` subcommand for quality control plots.

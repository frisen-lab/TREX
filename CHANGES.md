# TREX Changelog

## development version

* #7: Set default Jaccard threshold to 0.7 so that it matches R output.
* #33: Added a `trex qc` subcommand for quality control plots.
* #44: Fix some inconsistencies in the way cloneIDs are error-corrected.
* #30: Added doublet filtering. Cells that appear to be connected to two
  subclusters are detected and removed. The cell IDs of detected doublets
  is written to `doublets.txt`. Doublet detection can be disabled with
  `--keep-doublets`.
* [#40](https://github.com/frisen-lab/TREX/issues/40): Add a filtering option for excluding certain cloneIDs
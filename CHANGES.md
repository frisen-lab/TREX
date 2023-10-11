# TREX Changelog

## development version

* #7: Set default Jaccard threshold to 0.7 so that it matches R output.
* #33: Added a `trex qc` subcommand for quality control plots.
* #44: Fixed some inconsistencies in the way cloneIDs are error-corrected.
* #30: Added doublet filtering. Cells that appear to be connected to two
  subclusters are detected and removed. The cell IDs of detected doublets
  are written to `doublets.txt`. Doublet detection can be disabled with
  `--keep-doublets`.
* Added a filter to `run10x` that always removes low-complexity cloneIDs
  (those consisting of a single, repeated nucleotide such as `AAAA...`)

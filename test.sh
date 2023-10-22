#!/bin/bash
set -euo pipefail

# Ensure it is installed
samtools --version > /dev/null

pytest tests

set -x

trex smartseq2 --delete --chr EGFP-30N --output trex_smartseq2_run --start 2330 --end 2359 --read-matrix --readcount-threshold 1 tests/data/smartseq2_test.bam
diff -u <(sed 1d tests/expected_smartseq2/log.txt) <(sed 1d trex_smartseq2_run/log.txt)

trex smartseq3 --delete --output trex_smartseq3_run --umi-matrix -s 2330 -e 2359 tests/data/smartseq3_test.bam
diff -u <(sed 1d tests/expected_smartseq3/log.txt) <(sed 1d trex_smartseq3_run/log.txt)

trex qc trex_run

# Ensure --per-cell works
trex run10x --per-cell --delete -s 695 -e 724 -o trex_per_cell tests/data
diff -u <(sed 1d tests/expected_per_cell/log.txt) <(sed 1d trex_per_cell/log.txt)

# Ensure --filter-clone-id works
trex run10x --per-cell --delete -s 695 -e 724 -o trex_filter_clone_id --filter-cloneids tests/data/exclude_list.csv tests/data
diff -u <(sed 1d tests/expected_filter_clone_id/log.txt) <(sed 1d trex_filter_clone_id/log.txt)

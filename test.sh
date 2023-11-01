#!/bin/bash
set -euo pipefail

# Ensure it is installed
samtools --version > /dev/null

pytest tests

set -x


trex qc trex_run

# Ensure --per-cell works
trex run10x --per-cell --delete -s 695 -e 724 -o trex_per_cell tests/data
diff -u <(sed 1d tests/expected_per_cell/log.txt) <(sed 1d trex_per_cell/log.txt)

# Ensure --filter-clone-id works
trex run10x --per-cell --delete -s 695 -e 724 -o trex_filter_clone_id --filter-cloneids tests/data/exclude_list.csv tests/data
diff -u <(sed 1d tests/expected_filter_clone_id/log.txt) <(sed 1d trex_filter_clone_id/log.txt)

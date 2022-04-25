#!/bin/bash
set -euo pipefail
pytest tests
trex run10x --delete --loom --umi-matrix -s 695 -e 724 tests/data/
diff -ur -xdata.loom -xlog.txt tests/expected trex_run
diff -u <(sed 1d tests/expected/log.txt) <(sed 1d trex_run/log.txt)

trex smartseq2 --delete --chr EGFP-30N --output trex_smartseq2_run --start 2330 --end 2359 --read-matrix --readcount-threshold 1 tests/data/smartseq2_test.bam
diff -u <(sed 1d tests/expected_smartseq2/log.txt) <(sed 1d trex_smartseq2_run/log.txt)

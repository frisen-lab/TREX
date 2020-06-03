#!/bin/bash
set -euo pipefail pytest tests
braintrace run10x --delete --loom --umi-matrix -s 695 -e 724 tests/data/
diff -ur -xdata.loom -xlog.txt tests/expected braintrace_run
diff <(sed 1d tests/expected/log.txt) <(sed 1d braintrace_run/log.txt)
braintrace smartseq3 --delete --chr EGFP-30N --output braintrace_smartseq3_run --start 2330 --end 2359 --read-matrix --readcount-threshold 1 tests/data/smartseq3_test.bam
diff <(sed 1d tests/expected_smartseq3/log.txt) <(sed 1d braintrace_smartseq3_run/log.txt)

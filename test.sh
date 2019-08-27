#!/bin/bash
set -euo pipefail
pytest tests
braintrace --delete --loom --umi-matrix -s 695 -e 724 tests/data/
diff -ur -xdata.loom -xlog.txt tests/expected braintrace_run
diff <(sed 1d tests/expected/log.txt) <(sed 1d braintrace_run/log.txt)

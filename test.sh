#!/bin/bash
set -euo pipefail
pytest tests
trex run10x --delete --loom --umi-matrix -s 695 -e 724 tests/data/
diff -ur -xdata.loom -xlog.txt tests/expected trex_run
diff -u <(sed 1d tests/expected/log.txt) <(sed 1d trex_run/log.txt)

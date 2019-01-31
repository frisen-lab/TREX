#!/bin/bash
set -euo pipefail
braintrace --delete --loom --no-plot -s 695 -e 724 tests/data/outs
diff -ur -xdata.loom tests/expected lineage_run && echo "Test ok" || ( echo "Test fail"; exit 1)

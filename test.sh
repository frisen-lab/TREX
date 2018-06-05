#!/bin/bash
set -euo pipefail
test -e lineage_run && rm -r lineage_run
braintrace -s 695 -e 724 mini
diff -ur mini-expected lineage_run && echo "Test ok" || ( echo "Test fail"; exit 1)

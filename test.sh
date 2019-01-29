#!/bin/bash
set -euo pipefail
braintrace --delete --no-plot -s 695 -e 724 mini
diff -ur mini-expected lineage_run && echo "Test ok" || ( echo "Test fail"; exit 1)

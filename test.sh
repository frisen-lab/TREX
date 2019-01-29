#!/bin/bash
set -euo pipefail
braintrace --delete -s 695 -e 724 mini
diff -ur -x graph.pdf mini-expected lineage_run && echo "Test ok" || ( echo "Test fail"; exit 1)

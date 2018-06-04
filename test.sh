#!/bin/bash
set -euo pipefail
test -e lineage_run && rm -r lineage_run
python3 180425_lineage_loom.py -p mini
diff -ur mini-expected lineage_run && echo "Test ok" || ( echo "Test fail"; exit 1)

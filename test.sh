#!/bin/bash
set -euo pipefail
pytest tests
braintrace --delete --loom --umi-matrix -s 695 -e 724 tests/data/
diff -ur -xdata.loom tests/expected braintrace_run && echo "Test ok" || ( echo "Test fail"; exit 1)

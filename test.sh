#!/bin/bash
set -euo pipefail

# Ensure it is installed
samtools --version > /dev/null

pytest tests

set -x

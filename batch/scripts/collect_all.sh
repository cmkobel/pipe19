#!/bin/bash

# This script should be supersed by the r script (collect_batches.r), right now, it has problems with the types, and I don't have time to fix it.

echo "catting all batches together ..."
cat output/*/*_integrated.tsv > integrated.tsv
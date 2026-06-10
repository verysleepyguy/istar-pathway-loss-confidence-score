#!/bin/bash
set -e

prefix="data/demo/"
pathways="data/pathways-demo.tsv"

./download_demo.sh "${prefix}"
./download_checkpoints.sh
./run_pathway.sh "${prefix}" "${pathways}"

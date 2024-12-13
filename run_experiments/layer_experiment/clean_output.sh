#!/bin/bash

for output_dir in ./output_*/; do
    echo "Cleaning $output_dir"
    find "$output_dir" -type f ! -name '*.aln' -exec rm -v {} \;
    find "$output_dir" -type d -mindepth 1 -exec rm -rv {} \;
done
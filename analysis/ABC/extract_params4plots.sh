#!/bin/bash

# Extracts parameter columns from ABC simulation results
# Creates rsc_100k_params4PEplots_new.csv from rsc_100k_subset_new.csv
# 
# Usage: bash extract_params4plots.sh [input_file] [output_file]

# Set default file names
INPUT_FILE="${1:-rsc_100k_subset_new.csv}"
OUTPUT_FILE="${2:-rsc_100k_params4PEplots_new.csv}"


# Extract columns 8-28 (N.LA through T.Grp5) using cut

cut -d',' -f8-28 "$INPUT_FILE" > "$OUTPUT_FILE"

# Get row count for verification
TOTAL_ROWS=$(wc -l < "$OUTPUT_FILE")

echo "Successfully created $OUTPUT_FILE"
echo "Total rows: $TOTAL_ROWS (including header)"


echo ""
echo "File ready for use with abc_posteriors_plots.R"
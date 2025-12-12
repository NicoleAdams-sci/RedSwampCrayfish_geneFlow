#!/bin/bash

# Script to combine all rsc_100k_CV_*.csv files into a single file
# ABC cross-validation results
# Usage: bash ../SHELL/combine_cv_files.sh


# Set working directory and output file
WORK_DIR="/mnt/research/Scribner_Lab/projects/RedSwampCrayfish_MISGP/OUT/rscABC/OUT2"
OUTPUT_FILE="rsc_100k_CV_combo_new.csv"

cd "$WORK_DIR" || { echo "Error: Cannot change to directory $WORK_DIR"; exit 1; }

echo "Combining CV files in directory: $(pwd)"
echo "Output file: $OUTPUT_FILE"

# Remove output file if it exists
if [ -f "$OUTPUT_FILE" ]; then
    echo "Removing existing $OUTPUT_FILE"
    rm "$OUTPUT_FILE"
fi

# Initialize variables
HEADER_WRITTEN=false
TOTAL_FILES=0
TOTAL_ROWS=0

# Function to process files
process_files() {
    local pattern=$1
    local method=$2
    local tolerance=$3
    
    echo "Processing $pattern files..."
    
    for file in $pattern; do
        if [ -f "$file" ]; then
            TOTAL_FILES=$((TOTAL_FILES + 1))
            
            # Extract model name and node from filename
            # Pattern: rsc_100k_CV_[method]_[model_name]-[tolerance]-[node].csv
            if [[ "$file" =~ rsc_100k_CV_${method}_(.+)-${tolerance}-([0-9]+)\.csv ]]; then
                MODEL="${BASH_REMATCH[1]}"
                NODE="${BASH_REMATCH[2]}"
                echo "  Parsed: method=$method, model=$MODEL, tolerance=$tolerance, node=$NODE"
            else
                # Alternative parsing method - split on underscores and dashes
                filename=$(basename "$file" .csv)
                # Remove the prefix part: rsc_100k_CV_[method]_
                temp="${filename#rsc_100k_CV_${method}_}"
                # Extract node (last part after final dash)
                NODE="${temp##*-}"
                # Extract tolerance (second to last part)
                temp2="${temp%-*}"
                TOLERANCE_PART="${temp2##*-}"
                # Everything before tolerance is the model name
                MODEL="${temp2%-*}"
                
                echo "  Alternative parsing: model=$MODEL, tolerance=$TOLERANCE_PART, node=$NODE"
                
                # Validate the parsing worked
                if [[ "$TOLERANCE_PART" == "$tolerance" && "$NODE" =~ ^[0-9]+$ ]]; then
                    echo "  Successfully parsed with alternative method"
                else
                    echo "Warning: Could not parse filename $file"
                    MODEL="unknown"
                    NODE="unknown"
                fi
            fi
            
            # Process the file
            if [ "$HEADER_WRITTEN" = false ]; then
                # First file: write header with additional columns
                {
                    head -n 1 "$file" | sed 's/$/,method,tol,node,source_file/'
                } >> "$OUTPUT_FILE"
                HEADER_WRITTEN=true
                echo "Header written from: $file"
            fi
            
            # Write data rows (skip header) with metadata
            tail -n +2 "$file" | while IFS= read -r line; do
                echo "$line,$method,$tolerance,$NODE,$file" >> "$OUTPUT_FILE"
                TOTAL_ROWS=$((TOTAL_ROWS + 1))
            done
            
            # Count rows in this file
            FILE_ROWS=$(tail -n +2 "$file" | wc -l)
            echo "  $file: $FILE_ROWS rows (model: $MODEL, node: $NODE)"
            
        fi
    done
}

# Process neural network files (tolerance 0.025)
if ls rsc_100k_CV_neuralnet_*.csv 1> /dev/null 2>&1; then
    process_files "rsc_100k_CV_neuralnet_*.csv" "nnet" "0.025"
else
    echo "No neural network CV files found"
fi

# Process multinomial logistic files (tolerance 0.005)  
if ls rsc_100k_CV_mnlogistic_*.csv 1> /dev/null 2>&1; then
    process_files "rsc_100k_CV_mnlogistic_*.csv" "mnlog" "0.005"
else
    echo "No multinomial logistic CV files found"
fi

# Final summary
if [ -f "$OUTPUT_FILE" ]; then
    FINAL_ROWS=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
    echo ""
    echo "=== SUMMARY ==="
    echo "Files processed: $TOTAL_FILES"
    echo "Total data rows: $FINAL_ROWS"
    echo "Output file: $OUTPUT_FILE"
    echo "File size: $(du -h "$OUTPUT_FILE" | cut -f1)"
    
    # Show breakdown by method
    echo ""
    echo "Breakdown by method:"
    if command -v awk >/dev/null 2>&1; then
        tail -n +2 "$OUTPUT_FILE" | awk -F',' '{print $(NF-3)}' | sort | uniq -c
    fi
    
    echo ""
    echo "First few rows of combined file:"
    head -n 5 "$OUTPUT_FILE"
    
else
    echo "Error: No output file was created"
    exit 1
fi

echo ""
echo "Script completed successfully!"
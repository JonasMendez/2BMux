#!/bin/bash

# ##############################################################################
#
# merge_fastq_pairs.sh
#
# A script to find paired-end FASTQ files (R1/R2) in a directory and
# concatenate them into a single merged file.
#
# Usage:
#   2. Make sure it's executable using: chmod +x merge_fastq_pairs.sh
#   3. Run it in the directory containing your FASTQ files: bash concatenatepairs.sh
#
# ##############################################################################

# --- Configuration ---
# Use the current directory as the target folder.
# To specify a different folder, change this to: FOLDER_PATH="/path/to/your/fastq/files"
FOLDER_PATH="."

# --- Script Body ---
echo "Starting FASTQ merge process in: $FOLDER_PATH"
echo "-------------------------------------------------"

# Use a 'while' loop with 'find' to safely handle filenames with spaces or special characters.
# We look for all R1 files and then derive the R2 and output filenames from them.
find "$FOLDER_PATH" -type f -name "*_R1_*.fastq*" | while read -r r1_file; do

    # Derive R2 and merged filenames from the R1 file path
    # This replaces the first occurrence of "_R1_" with "_R2_" or "_merged_"
    r2_file="${r1_file/_R1_/_R2_}"
    
    # For the output file, we remove the _R1_ part completely.
    # This handles common naming conventions like 'sample_L001_R1_001.fastq.gz'
    # becoming 'sample_L001_001.fastq.gz'.
    merged_file="${r1_file/_R1_/}"

    # Check if the corresponding R2 file actually exists
    if [ -f "$r2_file" ]; then
        # Check if the merged file already exists to avoid re-doing work
        if [ -f "$merged_file" ]; then
            echo "SKIPPED: Merged file '$merged_file' already exists."
        else
            # Announce the merge operation
            echo "MERGING:"
            echo "  R1: $r1_file"
            echo "  R2: $r2_file"
            echo "  Output: $merged_file"
            
            # Perform the concatenation
            cat "$r1_file" "$r2_file" > "$merged_file"
            
            echo "  -> Merge complete."
            echo "" # Add a blank line for readability
        fi
    else
        # Report if an R1 file is found without a matching R2
        echo "WARNING: Found R1 file but no matching R2 file:"
        echo "  R1: $r1_file"
        echo "  Expected R2: $r2_file"
        echo ""
    fi
done

echo "-------------------------------------------------"
echo "All paired files have been processed."


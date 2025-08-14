#!/usr/bin/env python3

import argparse
import csv
import sys

def create_sample_sheet(plate_number, matrix_path, ids_path):
    """
    Generates a sample sheet by mapping sample IDs to their plate, row, and barcode.
    Handles variable column counts and ignores empty or invalid cells in the matrix.
    """
    # --- 1. Define the Barcode Mapping ---
    barcode_map = {
        1: 'ACAC', 2: 'GTCT', 3: 'TGGT', 4: 'CACT',
        5: 'GATG', 6: 'TCAC', 7: 'CTGA', 8: 'AAGC',
        9: 'GTAG', 10: 'GACA', 11: 'GTGA', 12: 'AGTC'
    }

    # --- 2. Load the SampleID-to-Index Mapping ---
    sample_id_map = {}
    try:
        with open(ids_path, mode='r', encoding='utf-8') as infile:
            reader = csv.DictReader(infile)
            if 'SampleID' not in reader.fieldnames or 'Index' not in reader.fieldnames:
                print(f"Error: The ID file '{ids_path}' must have 'SampleID' and 'Index' headers.", file=sys.stderr)
                sys.exit(1)
            for row in reader:
                sample_id_map[row['Index']] = row['SampleID']
    except FileNotFoundError:
        print(f"Error: The specified ID file was not found: {ids_path}", file=sys.stderr)
        sys.exit(1)

    # --- 3. Process the Plate Matrix and Generate Output Data ---
    output_data = []
    try:
        with open(matrix_path, mode='r', encoding='utf-8') as infile:
            matrix_reader = csv.reader(infile)
            for row_num, row_contents in enumerate(matrix_reader, 1):
                if row_num > 8:
                    print(f"Warning: Matrix file '{matrix_path}' has more than 8 rows. Ignoring extra rows.", file=sys.stderr)
                    break
                
                for col_num, sample_index in enumerate(row_contents, 1):
                    # This is the key logic for ignoring cells.
                    # It checks three things:
                    # 1. Is the cell empty or just whitespace? (sample_index.strip())
                    # 2. Is the index from the cell actually in our ID map? (sample_id_map.get(...))
                    # If either is false, we skip this cell and move to the next one.
                    # This robustly handles empty cells, "Blank", "N/A", etc.
                    sample_id = sample_id_map.get(sample_index.strip())
                    if not sample_id:
                        # This cell is either empty, 'Blank', or an index not in the ID file. Ignore it.
                        continue

                    if col_num > 12:
                        print(f"Warning: Row {row_num} in '{matrix_path}' has more than 12 columns. Ignoring extra columns.", file=sys.stderr)
                        break
                    
                    barcode = barcode_map.get(col_num)
                    
                    output_data.append({
                        'Plate': plate_number,
                        'Row': row_num,
                        '3illBC': barcode,
                        'SampleID': sample_id
                    })

    except FileNotFoundError:
        print(f"Error: The specified matrix file was not found: {matrix_path}", file=sys.stderr)
        sys.exit(1)

    # --- 4. Write the Output CSV File ---
    output_filename = f"Plate{plate_number}_IDs.csv"
    try:
        with open(output_filename, mode='w', newline='', encoding='utf-8') as outfile:
            fieldnames = ['Plate', 'Row', '3illBC', 'SampleID']
            writer = csv.DictWriter(outfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(output_data)
        
        print(f"Successfully created sample sheet: {output_filename}")
        print(f"Total valid samples processed: {len(output_data)}")

    except IOError as e:
        print(f"Error: Could not write to output file '{output_filename}'. Reason: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create a sample sheet by mapping sample IDs to plate layout and barcodes.",
        epilog="Example: python MakeSampleIds.py -plate 1 -matrix plate1matrix.csv -ids samples.csv"
    )
    parser.add_argument("-plate", required=True, type=int, help="An integer (e.g., 1-9) for the plate number.")
    parser.add_argument("-matrix", required=True, type=str, help="Path to the plate matrix CSV (no headers).")
    parser.add_argument("-ids", required=True, type=str, help="Path to the SampleID-to-Index CSV (with headers).")
    args = parser.parse_args()

    create_sample_sheet(
        plate_number=args.plate,
        matrix_path=args.matrix,
        ids_path=args.ids
    )

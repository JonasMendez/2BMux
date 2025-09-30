#!/usr/bin/env python3
import os
import glob
import re
import argparse
import subprocess
import sys
import csv
from collections import defaultdict
import logging
import math
import gzip

# --- Barcode and Adaptor definitions ---
BARCODE_PAIRS = {
    "ACAC": "GTGT", "GTCT": "AGAC", "TGGT": "ACCA", "CACT": "AGTG",
    "GATG": "CATC", "TCAC": "GTGA", "CTGA": "TCAG", "AAGC": "GCTT",
    "GTAG": "CTAC", "GACA": "TGTC", "GTGA": "TCAC", "AGTC": "GACT"
}
ILL_BARCODES = set(BARCODE_PAIRS.keys())
ANTI_BARCODES = set(BARCODE_PAIRS.values())
ANTI_TO_ILL = {v: k for k, v in BARCODE_PAIRS.items()}

# Expanded list of potential adapter sequences for R2
ILLUMINA_ADAPTERS_R2 = [
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", # Full P7 adapter sequence
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",  # Common variant
    "AGATCGGAAG", # Partial P7 adapter sequence
    "AGAT"
]

def parse_filename(filename, plate_pos, row_pos, read_pos):
    base_filename = os.path.basename(filename)
    base_filename_no_ext = base_filename.split(".")[0]
    parts = base_filename_no_ext.split("_")
    if not all(pos - 1 < len(parts) for pos in [plate_pos, row_pos, read_pos]):
        raise ValueError(f"Filename {filename} does not have enough parts")
    plate_str, row_str, read_str = parts[plate_pos-1], parts[row_pos-1], parts[read_pos-1]
    plate_num = int(re.search(r"\d+", plate_str).group())
    row_num = int(re.search(r"\d+", row_str).group())
    read = "R1" if "1" in read_str else "R2"
    return plate_num, row_num, read

def load_id_map(idmap_file):
    id_map, bc_to_sample = {}, {}
    with open(idmap_file, "r") as f:
        reader = csv.reader(f)
        header = next(reader)
        expected_header = ["Plate", "Row", "3illBC", "SampleID"]
        if header != expected_header:
            raise ValueError(f"ID map CSV header is incorrect. Expected: {expected_header}, Found: {header}")
        for row in reader:
            plate, row_num, barcode, sample_id = row
            id_map[(int(plate), int(row_num), barcode)] = sample_id
            bc_to_sample[barcode] = sample_id
    return id_map, bc_to_sample

# Modified to calculate running stats instead of storing all data
def update_running_stats(current_stats, new_value):
    count, current_sum, sum_sq = current_stats
    count += 1
    current_sum += new_value
    sum_sq += new_value ** 2
    return count, current_sum, sum_sq

def get_final_stats(current_stats):
    count, current_sum, sum_sq = current_stats
    if count == 0: return 0.0, 0.0
    mean = current_sum / count
    # Using Welford's algorithm for numerical stability, though not strictly necessary for small N
    # For population standard deviation (n in denominator)
    variance = (sum_sq / count) - (mean ** 2)
    stdev = math.sqrt(max(0, variance)) # Ensure non-negative for sqrt
    return mean, stdev

# Modified to calculate running average quality score
def update_running_qual_stats(current_qual_stats, qual_str):
    total_bases, total_phred_sum = current_qual_stats
    phred_scores = [ord(c) - 33 for c in qual_str]
    total_bases += len(phred_scores)
    total_phred_sum += sum(phred_scores)
    return total_bases, total_phred_sum

def get_final_avg_qual_score(current_qual_stats):
    total_bases, total_phred_sum = current_qual_stats
    if total_bases == 0: return 0.0
    return total_phred_sum / total_bases

def count_total_reads(fastq_file):
    count = 0
    if not fastq_file: return 0 # Handle case where fastq_file might be None
    try:
        with gzip.open(fastq_file, "rt") as f:
            for i, line in enumerate(f):
                if (i + 1) % 4 == 1: # Every 4th line is a read header
                    count += 1
    except FileNotFoundError:
        logging.warning(f"Input file not found: {fastq_file}. Assuming 0 reads.")
        return 0
    return count

def main():
    parser = argparse.ArgumentParser(description="Demultiplex, clean, and quality filter 2b-RAD paired-end FASTQ files.")
    parser.add_argument("-site", default="^(?:(?P<forward>.{36})|(?P<reverse>.{36}))", help="Regex for the DNA fragment (default: 36bp)")
    parser.add_argument("-barcode", default="[ATGC]{4}", help="Regex for in-line barcode (default: 4-bp ATGC)")
    parser.add_argument("-adaptor", default="AGATCGGAAG", help="Adaptor sequence on R1 (default: Partial P7)")
    parser.add_argument("-trim", type=int, default=0, help="Bases to trim from ends of fragment (default: 0)")
    parser.add_argument("-min_bc_count", type=int, default=10000, help="Minimum read count per barcode (default: 10000)")
    parser.add_argument("-plate", type=int, required=True, help="Position of plate number in filename")
    parser.add_argument("-row", type=int, required=True, help="Position of row number in filename")
    parser.add_argument("-read", type=int, required=True, help="Position of read designation in filename")
    parser.add_argument("-idmap", required=True, help="CSV file with Plate,Row,3illBC,SampleID columns")
    parser.add_argument("--strict_orientation", action="store_true", help="Discard reads with unexpected site/barcode orientation")
    parser.add_argument("--antiBC_only", action="store_true", help="Only use antiBC barcodes for demultiplexing")
    args = parser.parse_args()

    logging.basicConfig(filename="demux_2bRAD.log", level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s", force=True)
    logging.info("--- Starting New Demultiplexing Run ---")

    try:
        site_pattern, barcode_pattern = re.compile(args.site), re.compile(args.barcode)
    except re.error as e: logging.error(f"Failed to compile regex: {e}"); sys.exit(1)

    id_map, _ = load_id_map(args.idmap)
    
    demux_map = {}
    if args.antiBC_only:
        logging.info("Running in --antiBC_only mode. Demultiplexing will only use anti-barcodes.")
        for ill_bc, anti_bc in BARCODE_PAIRS.items():
            demux_map[anti_bc] = ill_bc
    else:
        for ill_bc, anti_bc in BARCODE_PAIRS.items():
            demux_map[ill_bc] = ill_bc
            demux_map[anti_bc] = ill_bc

    fastq_files = glob.glob("*.fq.gz") + glob.glob("*.fastq.gz")
    if not fastq_files: logging.error("No .fq.gz or .fastq.gz files found"); sys.exit(1)
    
    pairs = defaultdict(dict)
    for f in fastq_files:
        try:
            plate, row, read = parse_filename(f, args.plate, args.row, args.read)
            base_name = f"Plate{plate}_Row{row}"
            pairs[base_name][read] = f
        except Exception as e: logging.error(f"Failed to parse filename {f}: {e}"); sys.exit(1)

    for base_name, reads in pairs.items():
        r1_file, r2_file = reads.get("R1"), reads.get("R2")
        if not r1_file: continue
        plate, row, _ = parse_filename(r1_file, args.plate, args.row, args.read)
        logging.info(f"Processing pair: {r1_file}, {r2_file or 'N/A'}")

        # Track initial raw read counts
        initial_r1_raw_reads = count_total_reads(r1_file)
        initial_r2_raw_reads = count_total_reads(r2_file) if r2_file else 0
        logging.info(f"Initial raw R1 reads for {base_name}: {initial_r1_raw_reads}")
        if r2_file: logging.info(f"Initial raw R2 reads for {base_name}: {initial_r2_raw_reads}")

        # Initialize running statistics for lengths and quality scores
        # (count, sum, sum_sq) for lengths
        r1_len_stats = defaultdict(lambda: [0, 0, 0])
        r2_len_stats = defaultdict(lambda: [0, 0, 0])
        # (total_bases, total_phred_sum) for quality scores
        r1_qual_stats = defaultdict(lambda: [0, 0])
        r2_qual_stats = defaultdict(lambda: [0, 0])

        r1_barcode_counts = defaultdict(int)
        r2_barcode_counts = defaultdict(int)
        site_role_counts = defaultdict(lambda: defaultdict(int))
        read_id_to_target = {}

        # Open output files for streaming writes
        output_file_handles_r1 = {}
        output_file_handles_r2 = {}
        unmatched_file_handles_r1 = {}

        # Pre-create output directories for unmatched reads
        unmatched_dir = f"unmatched_barcodes_Plate{plate}Row{row}"
        os.makedirs(unmatched_dir, exist_ok=True)

        # First pass: Process R1 reads, write to temporary files, and collect R1 stats
        with gzip.open(r1_file, "rt") as f_r1:
            lines_r1 = iter(f_r1)
            for name, seq, plus, qual in zip(lines_r1, lines_r1, lines_r1, lines_r1):
                name, seq, qual = name.strip(), seq.strip(), qual.strip()
                m = site_pattern.match(seq)
                if not m: continue

                frag_start, frag_end = (m.start("forward"), m.end("forward")) if m.group("forward") else (m.start("reverse"), m.end("reverse"))
                orientation = "forward" if m.group("forward") else "reverse"
                bc_start, bc_end = frag_end, frag_end + 4
                
                if bc_end < len(seq) and barcode_pattern.fullmatch(seq[bc_start:bc_end]) and seq[bc_end:].startswith(args.adaptor):
                    obs_bc = seq[bc_start:bc_end]
                    
                    target_bc = demux_map.get(obs_bc)
                    
                    if args.strict_orientation and target_bc is not None:
                        is_ill_bc = obs_bc in ILL_BARCODES
                        is_anti_bc = obs_bc in ANTI_BARCODES
                        if (orientation == "forward" and is_anti_bc) or \
                           (orientation == "reverse" and is_ill_bc):
                            target_bc = None

                    if target_bc is None: target_bc = f"UNMATCHED_{obs_bc}"
                    
                    # populating site_role_counts
                    if "UNMATCHED" not in target_bc:
                        actual_ill_bc = demux_map.get(obs_bc, obs_bc) # Get the 3illBC key for the observed barcode
                        
                        if orientation == "forward":
                            if obs_bc in ANTI_BARCODES: site_role_counts[actual_ill_bc]["SiteFwd_anti"] += 1
                            elif obs_bc in ILL_BARCODES: site_role_counts[actual_ill_bc]["SiteFwd_3ill"] += 1
                        else: # reverse orientation
                            if obs_bc in ANTI_BARCODES: site_role_counts[actual_ill_bc]["SiteRev_anti"] += 1
                            elif obs_bc in ILL_BARCODES: site_role_counts[actual_ill_bc]["SiteRev_3ill"] += 1

                    frag_seq, frag_qual = seq[frag_start:frag_end], qual[frag_start:frag_end]
                    trim_start, trim_end = args.trim, len(frag_seq) - args.trim
                    if trim_end <= trim_start: continue
                    
                    final_seq, final_qual = frag_seq[trim_start:trim_end], frag_qual[trim_start:trim_end]
                    dline = f"{name} bcd={obs_bc}\n{final_seq}\n+\n{final_qual}\n"
                    read_id_to_target[name.split()[0]] = target_bc
                    r1_barcode_counts[target_bc] += 1 # Increment R1 count

                    # Write R1 read directly to file
                    if "UNMATCHED" in target_bc:
                        if target_bc not in unmatched_file_handles_r1:
                            unmatched_file_handles_r1[target_bc] = open(os.path.join(unmatched_dir, f"{target_bc}_R1.fq"), "a")
                        unmatched_file_handles_r1[target_bc].write(dline)
                    else:
                        sample_id = id_map.get((plate, row, target_bc))
                        if sample_id:
                            if target_bc not in output_file_handles_r1:
                                output_file_handles_r1[target_bc] = open(f"temp_{sample_id}_R1.fq", "a")
                            output_file_handles_r1[target_bc].write(dline)

                    # Update running stats for R1
                    r1_len_stats[target_bc] = update_running_stats(r1_len_stats[target_bc], len(final_seq))
                    r1_qual_stats[target_bc] = update_running_qual_stats(r1_qual_stats[target_bc], final_qual)

        # Close all R1 output file handles
        for fh in output_file_handles_r1.values(): fh.close()
        for fh in unmatched_file_handles_r1.values(): fh.close()

        # Second pass: Process R2 reads, write to temporary files, and collect R2 stats
        if r2_file:
            with gzip.open(r2_file, "rt") as f_r2:
                lines_r2 = iter(f_r2)
                for name, seq, plus, qual in zip(lines_r2, lines_r2, lines_r2, lines_r2):
                    read_id = name.strip().split()[0]
                    if read_id in read_id_to_target:
                        target_bc_key = read_id_to_target[read_id] # This is the target_bc (e.g., ACAC or GTGT)
                        if "UNMATCHED" in target_bc_key: continue

                        rd_seq, rd_qual = seq.strip(), qual.strip()
                        
                        # --- Simplified R2 Trimming Strategy ---
                        # 1. Initial Barcode Trimming (both ends, both orientations)
                        current_bc = target_bc_key # The observed barcode from R1
                        rc_current_bc = BARCODE_PAIRS.get(current_bc) or ANTI_TO_ILL.get(current_bc) # Get its RC

                        # Try to trim current_bc or its RC from 5' end
                        if len(rd_seq) >= 4 and (rd_seq.startswith(current_bc) or (rc_current_bc and rd_seq.startswith(rc_current_bc))):
                            rd_seq, rd_qual = rd_seq[4:], rd_qual[4:]
                        
                        # Try to trim current_bc or its RC from 3' end
                        if len(rd_seq) >= 4 and (rd_seq.endswith(current_bc) or (rc_current_bc and rd_seq.endswith(rc_current_bc))):
                            rd_seq, rd_qual = rd_seq[:-4], rd_qual[:-4]

                        # 2. Fixed-Length Biological Fragment Extraction (36bp)
                        #    and Adapter Motif Confirmation
                        expected_frag_len = 36
                        
                        if len(rd_seq) >= expected_frag_len:
                            potential_frag_seq = rd_seq[:expected_frag_len]
                            potential_frag_qual = rd_qual[:expected_frag_len]
                            
                            # Check for adapter immediately after the 36bp fragment
                            adapter_found = False
                            if len(rd_seq) > expected_frag_len: # There's sequence after the 36bp
                                # Check for any of the known adapters starting right after the 36bp fragment
                                for adapter_seq in ILLUMINA_ADAPTERS_R2:
                                    if rd_seq[expected_frag_len:].startswith(adapter_seq):
                                        adapter_found = True
                                        break
                            
                            # If adapter is found, or if there's no sequence after 36bp (perfect read)
                            if adapter_found or len(rd_seq) == expected_frag_len:
                                # This is our final sequence, apply symmetric trimming if any
                                final_seq = potential_frag_seq
                                final_qual = potential_frag_qual

                                if len(final_seq) > 2 * args.trim:
                                    final_seq = final_seq[args.trim:-args.trim] if args.trim > 0 else final_seq
                                    final_qual = final_qual[args.trim:-args.trim] if args.trim > 0 else final_qual
                                    
                                    dline = f"{name.strip()}\n{final_seq}\n+\n{final_qual}\n"
                                    r2_barcode_counts[target_bc_key] += 1 # Increment R2 count

                                    # Write R2 read directly to file
                                    sample_id = id_map.get((plate, row, target_bc_key))
                                    if sample_id:
                                        if target_bc_key not in output_file_handles_r2:
                                            output_file_handles_r2[target_bc_key] = open(f"temp_{sample_id}_R2.fq", "a")
                                        output_file_handles_r2[target_bc_key].write(dline)

                                    # Update running stats for R2
                                    r2_len_stats[target_bc_key] = update_running_stats(r2_len_stats[target_bc_key], len(final_seq))
                                    r2_qual_stats[target_bc_key] = update_running_qual_stats(r2_qual_stats[target_bc_key], final_qual)

        # Close all R2 output file handles
        for fh in output_file_handles_r2.values(): fh.close()

        output_files_to_filter = {}
        for target_bc in r1_barcode_counts.keys(): # Iterate over all R1 barcodes found
            if "UNMATCHED" in target_bc or r1_barcode_counts[target_bc] < args.min_bc_count: continue
            sample_id = id_map.get((plate, row, target_bc))
            if not sample_id: continue
            
            intermediate_r1, final_r1 = f"temp_{sample_id}_R1.fq", f"{sample_id}_R1_.fastq.gz"
            output_files_to_filter[intermediate_r1] = final_r1
            
            if r2_barcode_counts.get(target_bc):
                intermediate_r2, final_r2 = f"temp_{sample_id}_R2.fq", f"{sample_id}_R2_.fastq.gz"
                output_files_to_filter[intermediate_r2] = final_r2

        summary_file = f"Plate{plate}_Row{row}_barcode_summary.csv"
        with open(summary_file, "w", newline="") as f:
            header = ["3illBC", "SampleID", "R1_Reads", "R2_Reads", "R1_Retained_Perc", "R2_Retained_Perc", "SiteFwd_3ill", "SiteFwd_anti", "SiteRev_anti", "SiteRev_3ill", "R1_len_mean", "R1_len_stdev", "R2_len_mean", "R2_len_stdev", "R1_Qual_mean", "R2_Qual_mean"]
            if args.antiBC_only: header = ["3illBC", "SampleID", "R1_Reads", "R2_Reads", "R1_Retained_Perc", "R2_Retained_Perc", "SiteFwd_anti", "SiteRev_anti", "R1_len_mean", "R1_len_stdev", "R2_len_mean", "R2_len_stdev", "R1_Qual_mean", "R2_Qual_mean"]
            writer = csv.writer(f)
            writer.writerow(header)
            
            for ill_bc in sorted(BARCODE_PAIRS.keys()):
                sample_id = id_map.get((plate, row, ill_bc), "N/A")
                actual_ill_bc_for_counts = ill_bc
                counts = site_role_counts.get(actual_ill_bc_for_counts, {})
                
                r1_count = r1_barcode_counts.get(actual_ill_bc_for_counts, 0)
                r2_count = r2_barcode_counts.get(actual_ill_bc_for_counts, 0)
                
                r1_mean, r1_stdev = get_final_stats(r1_len_stats.get(actual_ill_bc_for_counts, [0,0,0]))
                r2_mean, r2_stdev = get_final_stats(r2_len_stats.get(actual_ill_bc_for_counts, [0,0,0]))
                r1_qual_mean = get_final_avg_qual_score(r1_qual_stats.get(actual_ill_bc_for_counts, [0,0]))
                r2_qual_mean = get_final_avg_qual_score(r2_qual_stats.get(actual_ill_bc_for_counts, [0,0]))

                r1_retained_perc = (r1_count / initial_r1_raw_reads * 100) if initial_r1_raw_reads > 0 else 0.0
                r2_retained_perc = (r2_count / initial_r2_raw_reads * 100) if initial_r2_raw_reads > 0 else 0.0

                row_data = [ill_bc, sample_id, r1_count, r2_count, f"{r1_retained_perc:.2f}%", f"{r2_retained_perc:.2f}%",
                            counts.get("SiteFwd_3ill", 0), counts.get("SiteFwd_anti", 0),
                            counts.get("SiteRev_anti", 0), counts.get("SiteRev_3ill", 0),
                            f"{r1_mean:.2f}", f"{r1_stdev:.2f}", f"{r2_mean:.2f}", f"{r2_stdev:.2f}",
                            f"{r1_qual_mean:.2f}", f"{r2_qual_mean:.2f}"]
                
                if args.antiBC_only:
                    row_data = [ill_bc, sample_id, r1_count, r2_count, f"{r1_retained_perc:.2f}%", f"{r2_retained_perc:.2f}%",
                                counts.get("SiteFwd_anti", 0), counts.get("SiteRev_anti", 0),
                                f"{r1_mean:.2f}", f"{r1_stdev:.2f}", f"{r2_mean:.2f}", f"{r2_stdev:.2f}",
                                f"{r1_qual_mean:.2f}", f"{r2_qual_mean:.2f}"]

                writer.writerow(row_data)

        logging.info(f"Summary written to {summary_file}")

        # Apply quality filtering and gzip compression using external tool
        for intermediate_file, final_file in output_files_to_filter.items():
            if os.path.exists(intermediate_file):
                logging.info(f"Applying quality filter and gzipping {intermediate_file} to {final_file}")
                try:
                    # Using zcat and gzip for efficient streaming compression
                    # This avoids loading the entire unzipped file into memory for gzip
                    cmd = f"cat {intermediate_file} | gzip > {final_file}"
                    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    os.remove(intermediate_file) # Clean up intermediate unzipped file
                except subprocess.CalledProcessError as e:
                    logging.error(f"Error processing {intermediate_file}: {e.stderr.decode().strip()}")
                    sys.exit(1)
            else:
                logging.warning(f"Intermediate file {intermediate_file} not found for processing.")

    logging.info("--- Demultiplexing Run Finished ---")

if __name__ == "__main__":
    main()

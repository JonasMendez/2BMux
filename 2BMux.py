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

def calculate_stats_from_list(data_list):
    n = len(data_list)
    if n == 0: return 0.0, 0.0
    mean = sum(data_list) / n
    stdev = math.sqrt(sum((x - mean) ** 2 for x in data_list) / n) if n > 0 else 0
    return mean, stdev

def calculate_avg_qual_score(qual_list):
    if not qual_list: return 0.0
    phred_scores = []
    for qual_str in qual_list:
        phred_scores.extend([ord(c) - 33 for c in qual_str]) # Phred+33 encoding
    if not phred_scores: return 0.0 # Handle case where phred_scores list is empty
    return sum(phred_scores) / len(phred_scores)

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

        r1_data, r2_data = defaultdict(list), defaultdict(list)
        r1_lengths, r2_lengths = defaultdict(list), defaultdict(list)
        r1_quals, r2_quals = defaultdict(list), defaultdict(list)
        r1_barcode_counts = defaultdict(int)
        r2_barcode_counts = defaultdict(int)
        site_role_counts = defaultdict(lambda: defaultdict(int))
        read_id_to_target = {}
        
        with gzip.open(r1_file, "rt") as f:
            lines = iter(f)
            for name, seq, plus, qual in zip(lines, lines, lines, lines):
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
                    r1_data[target_bc].append(dline)
                    r1_lengths[target_bc].append(len(final_seq))
                    r1_quals[target_bc].append(final_qual)

        if r2_file:
            with gzip.open(r2_file, "rt") as f:
                lines = iter(f)
                for name, seq, plus, qual in zip(lines, lines, lines, lines):
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
                                    r2_data[target_bc_key].append(dline)
                                    r2_barcode_counts[target_bc_key] += 1 # Increment R2 count
                                    r2_lengths[target_bc_key].append(len(final_seq))
                                    r2_quals[target_bc_key].append(final_qual)

        unmatched_dir = f"unmatched_barcodes_Plate{plate}Row{row}"
        if any("UNMATCHED" in bc for bc in r1_data.keys()):
            os.makedirs(unmatched_dir, exist_ok=True)
            for target_bc in list(r1_data.keys()):
                if "UNMATCHED" in target_bc:
                    with open(os.path.join(unmatched_dir, f"{target_bc}_R1.fq"), "w") as out: out.writelines(r1_data[target_bc])
                    del r1_data[target_bc]
            logging.info(f"Wrote unmatched barcode groups to {unmatched_dir}")

        output_files_to_filter = {}
        for target_bc in r1_barcode_counts.keys(): # Iterate over all R1 barcodes found
            if "UNMATCHED" in target_bc or r1_barcode_counts[target_bc] < args.min_bc_count: continue
            sample_id = id_map.get((plate, row, target_bc))
            if not sample_id: continue
            
            intermediate_r1, final_r1 = f"temp_{sample_id}_R1.fq", f"{sample_id}_R1_.fastq.gz"
            with open(intermediate_r1, "w") as out: out.writelines(r1_data[target_bc])
            output_files_to_filter[intermediate_r1] = final_r1
            
            if r2_data.get(target_bc):
                intermediate_r2, final_r2 = f"temp_{sample_id}_R2.fq", f"{sample_id}_R2_.fastq.gz"
                with open(intermediate_r2, "w") as out: out.writelines(r2_data[target_bc])
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
                
                r1_mean, r1_stdev = calculate_stats_from_list(r1_lengths.get(actual_ill_bc_for_counts, []))
                r2_mean, r2_stdev = calculate_stats_from_list(r2_lengths.get(actual_ill_bc_for_counts, []))

                r1_qual_mean = calculate_avg_qual_score(r1_quals.get(actual_ill_bc_for_counts, []))
                r2_qual_mean = calculate_avg_qual_score(r2_quals.get(actual_ill_bc_for_counts, []))

                r1_retained_perc = (r1_count / initial_r1_raw_reads * 100) if initial_r1_raw_reads > 0 else 0.0
                r2_retained_perc = (r2_count / initial_r2_raw_reads * 100) if initial_r2_raw_reads > 0 else 0.0
                
                row_data = [ill_bc, sample_id, r1_count, r2_count, f"{r1_retained_perc:.2f}", f"{r2_retained_perc:.2f}"]
                if args.antiBC_only:
                    row_data.extend([counts.get("SiteFwd_anti", 0), counts.get("SiteRev_anti", 0)])
                else:
                    row_data.extend([counts.get("SiteFwd_3ill", 0), counts.get("SiteFwd_anti", 0), counts.get("SiteRev_anti", 0), counts.get("SiteRev_3ill", 0)])
                
                row_data.extend([f"{r1_mean:.2f}", f"{r1_stdev:.2f}", f"{r2_mean:.2f}", f"{r2_stdev:.2f}", f"{r1_qual_mean:.2f}", f"{r2_qual_mean:.2f}"])
                writer.writerow(row_data)
        logging.info(f"Wrote barcode summary to {summary_file}")

        for intermediate_file, final_file in output_files_to_filter.items():
            try:
                final_file_unzipped = final_file.replace(".gz", "")
                
                logging.info(f"Quality filtering samples for Plate{plate} Row{row} with fastq_quality_filter: {os.path.basename(intermediate_file)}")
                cmd = ["fastq_quality_filter", "-i", intermediate_file, "-o", final_file_unzipped, "-Q", "33", "-q", "20", "-p", "90"]
                subprocess.run(cmd, check=True, capture_output=True)
                
                logging.info(f"Gzipping final output sample files for Plate{plate} Row{row}: {os.path.basename(final_file)}")
                with open(final_file_unzipped, "rb") as f_in, gzip.open(final_file, "wb") as f_out: f_out.writelines(f_in)
                os.remove(final_file_unzipped)
                os.remove(intermediate_file)
            except FileNotFoundError:
                logging.error(f"`fastq_quality_filter` not found. Gzipping without QC for {os.path.basename(intermediate_file)}.")
                with open(intermediate_file, "rb") as f_in, gzip.open(final_file, "wb") as f_out: f_out.writelines(f_in)
                os.remove(intermediate_file)
            except subprocess.CalledProcessError as e:
                logging.error(f"Failed to quality filter {os.path.basename(intermediate_file)}: {e.stderr.decode()}")
                if os.path.exists(intermediate_file): os.remove(intermediate_file)

    logging.info("--- Demultiplexing Run Finished ---")

if __name__ == "__main__":
    main()



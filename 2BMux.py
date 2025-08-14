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

# Hardcoded 3illBC and antiBC barcodes, these may be re hard-coded to your own specific barcodes
BARCODE_PAIRS = {
    "ACAC": "GTGT",
    "GTCT": "AGAC",
    "TGGT": "ACCA",
    "CACT": "AGTG",
    "GATG": "CATC",
    "TCAC": "GTGA",
    "CTGA": "TCAG",
    "AAGC": "GCTT",
    "GTAG": "CTAC",
    "GACA": "TGTC",
    "GTGA": "TCAC",
    "AGTC": "GACT"
}
ILL_BARCODES = set(BARCODE_PAIRS.keys())
ANTI_BARCODES = set(BARCODE_PAIRS.values())
ANTI_TO_ILL = {v: k for k, v in BARCODE_PAIRS.items()}

# Common Illumina adaptor sequences (partial); you may modify and add hardcoded adaptor sequences if needed.
ADAPTORS = [
    'AATGATACGGCG',  # Partial P5 adaptor
    'AGATCGGAAG',    # Partial P7 adaptor
    'CAAGCAGAAGAC'   # Partial P7 adaptor
]


def parse_filename(filename, plate_pos, row_pos, read_pos):
    #Parse plate number, row number, and read designation from filename
    # Extract filename without path and extension
    base_filename = os.path.basename(filename)
    base_filename_no_ext = base_filename.split('.')[0] # Remove .fastq.gz or .fq.gz

    parts = base_filename_no_ext.split('_')
    if plate_pos - 1 >= len(parts) or row_pos - 1 >= len(parts) or read_pos - 1 >= len(parts):
        raise ValueError(f"Filename {filename} does not have enough underscore-separated parts")

    # Plate number: Extract first integer > 0
    plate_str = parts[plate_pos - 1]
    plate_nums = re.findall(r'\d+', plate_str)
    plate_nums = [int(n) for n in plate_nums if int(n) > 0]
    if len(plate_nums) != 1:
        raise ValueError(f"Plate ID {plate_str} must contain exactly one integer > 0")

    # Row number: Extract first integer > 0
    row_str = parts[row_pos - 1]
    row_nums = re.findall(r'\d+', row_str)
    row_nums = [int(n) for n in row_nums if int(n) > 0]
    if len(row_nums) != 1:
        raise ValueError(f"Row ID {row_str} must contain exactly one integer > 0")

    # Read designation: Look for 1 or 2
    read_str = parts[read_pos - 1]
    read_match = re.search(r'[12]', read_str)
    if not read_match:
        raise ValueError(f"Read designation {read_str} must contain '1' or '2'")
    read = 'R1' if read_match.group(0) == '1' else 'R2'

    return plate_nums[0], row_nums[0], read


def load_id_map(idmap_file):
    #Load ID mapping from CSV
    id_map = {}
    bc_to_sample = {}
    with open(idmap_file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        if header != ['Plate', 'Row', '3illBC', 'SampleID']:
            raise ValueError("ID map CSV must have columns: Plate,Row,3illBC,SampleID")
        for row in reader:
            if len(row) != 4:
                raise ValueError(f"Invalid ID map row: {row}")
            plate, row_num, barcode, sample_id = row
            plate = int(plate)
            row_num = int(row_num)
            if barcode not in BARCODE_PAIRS:
                raise ValueError(f"3illBC {barcode} in ID map not in real barcodes")
            id_map[(plate, row_num, barcode)] = sample_id
            bc_to_sample[barcode] = sample_id
    return id_map, bc_to_sample


def main():
    parser = argparse.ArgumentParser(description="Demultiplex and quality filter 2b-RAD paired-end FASTQ files.")
    parser.add_argument(
        '-site',
        default='^(?:(?P<forward>.{12}CGA.{6}TGC.{12})|(?P<reverse>.{12}GCA.{6}TCG.{12}))',
        help='Regex pattern for restriction site (default: BcgI 2b-RAD)'
    )
    parser.add_argument('-barcode', default='[ATGC]{4}',
                        help='Regex pattern for in-line barcode (default: 4-bp ATGC)')
    parser.add_argument('-adaptor', default='|'.join(ADAPTORS),
                        help='Adaptor sequences used to trim R2 (default: Illumina P5/P7 partial sequences)')
    parser.add_argument('-trim', type=int, default=0,
                        help='Number of bases to trim from ends of fragment (default: 0)')
    parser.add_argument('-min_bc_count', type=int, default=10000,
                        help='Minimum read count per barcode for output (default: 10000)')
    parser.add_argument('-plate', type=int, required=True,
                        help='Position of plate number in filename (1-based)')
    parser.add_argument('-row', type=int, required=True,
                        help='Position of row number in filename (1-based)')
    parser.add_argument('-read', type=int, required=True,
                        help='Position of read designation (1 or 2) in filename (1-based)')
    parser.add_argument('-idmap', required=True,
                        help='CSV file with Plate,Row,3illBC,SampleID columns')
    parser.add_argument('--strict_orientation', action='store_true',
                        help='If set, discard reads where site orientation and barcode role disagree (default: allow but flag)')
    parser.add_argument('--antiBC_only', action='store_true',
                        help='If set, only consider antiBC barcodes for demultiplexing, regardless of site orientation.')
    parser.add_argument('--allR1', action='store_true',
                        help='If set, keep demultiplexed R1 reads even if they do not have an R2 pair.')

    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(
        filename='demux_2bRAD.log',
        level=logging.DEBUG,
        format='%(asctime)s - %(levelname)s - %(message)s',
        force=True
    )
    logging.debug("Logging initialized")

    # Check for barcode overlaps
    overlaps = [bc for bc in BARCODE_PAIRS if bc in ANTI_BARCODES]
    if overlaps:
        logging.warning(f"Barcode overlaps detected (3illBC present as antiBC too): {sorted(overlaps)}. Orientation must be used to disambiguate.")

    logging.info("Starting demultiplexing and quality filtering")

    # Compile regex patterns
    try:
        site_pattern = re.compile(args.site)
        barcode_pattern = re.compile(args.barcode)
        adaptor_pattern = re.compile(f'^(.*)(?:{args.adaptor}).*$')  # use CLI-provided adaptor pattern for R2 trimming
    except re.error as e:
        logging.error(f"Failed to compile regex patterns: {e}")
        print(f"Error: Failed to compile regex patterns: {e}", file=sys.stderr)
        sys.exit(1)

    # Load ID mapping
    try:
        id_map, bc_to_sample = load_id_map(args.idmap)
    except Exception as e:
        logging.error(f"Failed to load ID map: {e}")
        print(f"Error: Failed to load ID map: {e}", file=sys.stderr)
        sys.exit(1)

    # Find FASTQ files
    fastq_files = glob.glob('*.fq.gz') + glob.glob('*.fastq.gz')
    if not fastq_files:
        logging.error("No .fq.gz or .fastq.gz files found in current directory")
        print("Error: No .fq.gz or .fastq.gz files found in current directory", file=sys.stderr)
        sys.exit(1)

    # Pair R1 and R2 files
    pairs = {}
    for f in fastq_files:
        try:
            plate, row, read = parse_filename(f, args.plate, args.row, args.read)
            # Construct a unique base_name for pairing, excluding the read designation part
            # This assumes the part before the read designation is consistent for R1/R2 pairs
            # Example: plate1_row1_R1.fastq.gz -> plate1_row1.fastq.gz
            #          plate1_row1_R2.fastq.gz -> plate1_row1.fastq.gz
            parts_no_ext = os.path.basename(f).split('.')[0].split('_')
            base_name_parts = parts_no_ext[:args.read-1] + parts_no_ext[args.read:]
            base_name = '_'.join(base_name_parts)

            pairs.setdefault(base_name, {})[read] = f
            logging.debug(f"Parsed filename {f}: Plate={plate}, Row={row}, Read={read}")
        except Exception as e:
            logging.error(f"Failed to parse filename {f}: {e}")
            print(f"Error: Failed to parse filename {f}: {e}", file=sys.stderr)
            sys.exit(1)

    for base_name, reads in pairs.items():
        if 'R1' not in reads or ('R2' not in reads and not args.allR1):
            if 'R1' not in reads:
                logging.error(f"Missing R1 file for {base_name}")
                print(f"Error: Missing R1 file for {base_name}", file=sys.stderr)
            elif 'R2' not in reads and not args.allR1:
                logging.error(f"Missing R2 file for {base_name} and --allR1 is not set")
                print(f"Error: Missing R2 file for {base_name} and --allR1 is not set", file=sys.stderr)
            sys.exit(1)

    # Process each pair
    for base_name, reads in pairs.items():
        r1_file = reads['R1']
        r2_file = reads.get('R2') # R2 might be None if --allR1 is set and no R2 exists
        try:
            plate, row, _ = parse_filename(r1_file, args.plate, args.row, args.read)
        except Exception as e:
            logging.error(f"Failed to parse filename {r1_file}: {e}")
            print(f"Error: Failed to parse filename {r1_file}: {e}", file=sys.stderr)
            sys.exit(1)

        logging.info(f"Processing pair: {r1_file}, {r2_file if r2_file else 'N/A'}")

        # Demultiplex R1 and store barcode assignments
        r1_data = defaultdict(list)
        r2_data = defaultdict(list)
        barcode_counts = defaultdict(int)  # per target 3illBC
        # Detailed counts
        site_role_counts = defaultdict(lambda: defaultdict(int))  # target_3illBC -> {'site_fwd_3ill': n, 'site_fwd_anti': n, 'site_rev_anti': n, 'site_rev_3ill': n}

        read_id_to_target = {}  # read id -> target 3illBC
        r2_skipped = 0
        r2_no_adaptor = 0
        unmatched_r1_count = 0

        # Process R1 file
        try:
            process = subprocess.Popen(['gzip', '-dc', r1_file], stdout=subprocess.PIPE, text=True)
            lines = iter(process.stdout)
            while True:
                try:
                    name = next(lines).strip()
                    seq = next(lines).strip()
                    plus = next(lines).strip()
                    qual = next(lines).strip()
                    if not name.startswith('@') or plus != '+':
                        raise ValueError("Invalid FASTQ format")
                    m = site_pattern.match(seq)
                    # Enforce barcode immediately after site and adaptor after barcode on R1
                    if m:
                        # Exactly one of these will be non-None
                        frag = m.group('forward') if m.group('forward') else m.group('reverse')
                        orientation = 'forward' if m.group('forward') else 'reverse'
                        bc_start = m.end()
                        bc_end = bc_start + 4
                        if bc_end <= len(seq) and barcode_pattern.fullmatch(seq[bc_start:bc_end]) and seq[bc_end:].startswith('AGATCGGAAG'):
                            obs_bc = seq[bc_start:bc_end]
                            
                            target_bc = None
                            if args.antiBC_only:
                                # AntiBC-only mode: only assign if observed barcode is an antiBC
                                if obs_bc in ANTI_BARCODES:
                                    target_bc = ANTI_TO_ILL[obs_bc]
                                    # Log for summary, but all reads here are treated as 'antiBC' role
                                    if orientation == 'forward':
                                        site_role_counts[target_bc]['site_fwd_anti'] += 1
                                    else: # orientation == 'reverse'
                                        site_role_counts[target_bc]['site_rev_anti'] += 1
                                else:
                                    # In antiBC-only mode, if it's not an antiBC, it's unmatched
                                    target_bc = None # Will be handled as unmatched below
                            else:
                                # Original logic: Decide mapping based on site orientation and barcode role
                                if orientation == 'forward':
                                    # Expected: observed bc is 3illBC
                                    if obs_bc in ILL_BARCODES:
                                        target_bc = obs_bc
                                        site_role_counts[target_bc]['site_fwd_3ill'] += 1
                                    elif obs_bc in ANTI_BARCODES:
                                        target_bc = ANTI_TO_ILL[obs_bc]
                                        if args.strict_orientation:
                                            # discard unexpected combo
                                            continue
                                        site_role_counts[target_bc]['site_fwd_anti'] += 1
                                    else:
                                        # Unknown barcode
                                        target_bc = None
                                    
                                else:  # orientation == 'reverse'
                                    # Expected: observed bc is antiBC
                                    if obs_bc in ANTI_BARCODES:
                                        target_bc = ANTI_TO_ILL[obs_bc]
                                        site_role_counts[target_bc]['site_rev_anti'] += 1
                                    elif obs_bc in ILL_BARCODES:
                                        target_bc = obs_bc
                                        if args.strict_orientation:
                                            continue
                                        site_role_counts[target_bc]['site_rev_3ill'] += 1
                                    else:
                                        target_bc = None
                            
                            if target_bc is None:
                                # Record unmatched into special bucket (write later)
                                # Keep only R1 unmatched for inspection
                                trimmed_len = len(frag) - 2*args.trim if len(frag) >= 2*args.trim else 0
                                if trimmed_len > 0:
                                    rd_seq = frag[args.trim: len(frag) - args.trim]
                                    rd_qual = qual[args.trim: len(frag) - args.trim]
                                    dline = f"{name} bcd={obs_bc}\n{rd_seq}\n+\n{rd_qual}\n"
                                    r1_data[obs_bc].append(dline)
                                    unmatched_r1_count += 1
                                continue

                            # Build demultiplexed R1 record
                            trimmed_len = len(frag) - 2*args.trim if len(frag) >= 2*args.trim else 0
                            if trimmed_len <= 0:
                                continue
                            rd_seq = frag[args.trim: len(frag) - args.trim]
                            rd_qual = qual[args.trim: len(frag) - args.trim]
                            dline = f"{name}\n{rd_seq}\n+\n{rd_qual}\n"
                            read_id = name.split()[0]
                            read_id_to_target[read_id] = target_bc
                            barcode_counts[target_bc] += 1
                            r1_data[target_bc].append(dline)
                        else:
                            # No valid 4bp barcode + adaptor
                            continue
                    else:
                        # No site match
                        continue
                except StopIteration:
                    break
            process.stdout.close()
            process.wait()
            if process.returncode != 0:
                raise subprocess.CalledProcessError(process.returncode, 'gzip -dc')
        except Exception as e:
            logging.error(f"Failed to process R1 file {r1_file}: {e}")
            print(f"Error: Failed to process R1 file {r1_file}: {e}", file=sys.stderr)
            sys.exit(1)

        # Process R2 file (only if r2_file exists)
        if r2_file:
            try:
                process = subprocess.Popen(['gzip', '-dc', r2_file], stdout=subprocess.PIPE, text=True)
                lines = iter(process.stdout)
                while True:
                    try:
                        name = next(lines).strip()
                        seq = next(lines).strip()
                        plus = next(lines).strip()
                        qual = next(lines).strip()
                        if not name.startswith('@') or plus != '+':
                            raise ValueError("Invalid FASTQ format")
                        read_id = name.split()[0]
                        if read_id in read_id_to_target:
                            # Check for adaptor sequences for trimming R2
                            m = adaptor_pattern.match(seq)
                            if m:
                                rd_seq = m.group(1)
                                rd_qual = qual[:len(rd_seq)]
                            else:
                                rd_seq = seq
                                rd_qual = qual
                                r2_no_adaptor += 1
                            # Apply trimming
                            if len(rd_seq) >= 2 * args.trim:
                                rd_seq = rd_seq[args.trim: len(rd_seq) - args.trim]
                                rd_qual = rd_qual[args.trim: len(rd_qual) - args.trim]
                            else:
                                r2_skipped += 1
                                continue
                            dline = f"{name}\n{rd_seq}\n+\n{rd_qual}\n"
                            target_bc = read_id_to_target[read_id]
                            r2_data[target_bc].append(dline)
                        else:
                            r2_skipped += 1
                            continue
                    except StopIteration:
                        break
                process.stdout.close()
                process.wait()
                if process.returncode != 0:
                    raise subprocess.CalledProcessError(process.returncode, 'gzip -dc')
            except Exception as e:
                logging.error(f"Failed to process R2 file {r2_file}: {e}")
                print(f"Error: Failed to process R2 file {r2_file}: {e}", file=sys.stderr)
                sys.exit(1)
            logging.info(f"Skipped {r2_skipped} R2 reads due to missing R1 pair or insufficient length")
            logging.info(f"Processed {r2_no_adaptor} R2 reads without detected adaptor sequences")

        # Write demultiplexed files
        unmatched_dir = f'unmatched_barcodes_Plate{plate}Row{row}'
        os.makedirs(unmatched_dir, exist_ok=True)
        unmatched_written_count = 0
        for target_bc in list(r1_data.keys()):
            # Unknown barcodes are not in ILL/ANTI; keep for inspection only
            if target_bc not in ILL_BARCODES and target_bc not in ANTI_BARCODES:
                outname = os.path.join(unmatched_dir, f"Plate{plate}_Row{row}_{target_bc}_R1.fq")
                with open(outname, 'w') as out:
                    for dline in r1_data[target_bc]:
                        out.write(dline)
                outname_r2 = os.path.join(unmatched_dir, f"Plate{plate}_Row{row}_{target_bc}_R2.fq")
                with open(outname_r2, 'w') as out:
                    for dline in r2_data.get(target_bc, []):
                        out.write(dline)
                unmatched_written_count += 1
                # Remove to avoid further processing
                del r1_data[target_bc]
                if target_bc in r2_data:
                    del r2_data[target_bc]
        if unmatched_written_count > 0:
            logging.info(f"Wrote {unmatched_written_count} unmatched barcode groups to {unmatched_dir}")

        output_files = set()
        for target_bc in sorted(r1_data.keys()):
            if barcode_counts[target_bc] < args.min_bc_count:
                logging.info(f"Skipping barcode {target_bc} with {barcode_counts[target_bc]} reads (below min_bc_count {args.min_bc_count})")
                continue
            
            # Final output filename format: SampleID_R1_.fastq
            sample_id = id_map.get((plate, row, target_bc))
            if not sample_id:
                logging.warning(f"No sample ID found for Plate{plate}_Row{row}_{target_bc}. Skipping output for this barcode.")
                continue

            final_r1_out = f"{sample_id}_R1_.fastq"
            final_r2_out = f"{sample_id}_R2_.fastq"

            # Write intermediate files for quality filtering
            intermediate_r1_out = f"Plate{plate}_Row{row}_{target_bc}_R1.fq"
            with open(intermediate_r1_out, 'w') as out:
                for dline in r1_data[target_bc]:
                    out.write(dline)
            output_files.add(intermediate_r1_out)

            if r2_file:  # Only write R2 if an actual R2 file exists
                intermediate_r2_out = f"Plate{plate}_Row{row}_{target_bc}_R2.fq"
                with open(intermediate_r2_out, 'w') as out:
                    for dline in r2_data.get(target_bc, []):
                        out.write(dline)
                output_files.add(intermediate_r2_out)

            # Logging for intermediate files is now combined with quality filter step

        # Write summary table (per plate-row)
        summary_file = f"Plate{plate}_Row{row}_barcode_summary.csv"
        with open(summary_file, 'w', newline='') as f:
            writer = csv.writer(f)
            # Adjust header based on antiBC_only mode
            if args.antiBC_only:
                writer.writerow(['3illBC', 'SampleID', 'TotalReads', 'SiteFwd_anti', 'SiteRev_anti'])
            else:
                writer.writerow(['3illBC', 'SampleID', 'TotalReads', 'SiteFwd_3ill', 'SiteFwd_anti', 'SiteRev_anti', 'SiteRev_3ill'])
            
            for ill_bc, anti_bc in sorted(BARCODE_PAIRS.items()):
                sample_id = id_map.get((plate, row, ill_bc), bc_to_sample.get(ill_bc, 'Unknown'))
                counts = site_role_counts.get(ill_bc, {})
                
                if args.antiBC_only:
                    sf_anti = counts.get('site_fwd_anti', 0)
                    sr_anti = counts.get('site_rev_anti', 0)
                    total = sf_anti + sr_anti
                    writer.writerow([ill_bc, sample_id, total, sf_anti, sr_anti])
                else:
                    sf3 = counts.get('site_fwd_3ill', 0)
                    sfa = counts.get('site_fwd_anti', 0)
                    sra = counts.get('site_rev_anti', 0)
                    sr3 = counts.get('site_rev_3ill', 0)
                    total = sf3 + sfa + sra + sr3
                    writer.writerow([ill_bc, sample_id, total, sf3, sfa, sra, sr3])
        logging.info(f"Wrote barcode summary to {summary_file}")

        # Quality filter and rename
        for intermediate_outname in list(output_files): # Iterate over a copy as we modify the set
            try:
                # Determine if it's R1 or R2 and get corresponding final name
                if '_R1.fq' in intermediate_outname:
                    barcode = intermediate_outname.split('_')[-2]
                    sample_id = id_map.get((plate, row, barcode))
                    final_out = f"{sample_id}_R1_.fastq"
                elif '_R2.fq' in intermediate_outname:
                    barcode = intermediate_outname.split('_')[-2]
                    sample_id = id_map.get((plate, row, barcode))
                    final_out = f"{sample_id}_R2_.fastq"
                else:
                    logging.warning(f"Unexpected intermediate file name format: {intermediate_outname}. Skipping quality filter.")
                    continue

                cmd = [
                    'fastq_quality_filter',
                    '-i', intermediate_outname,
                    '-o', final_out,
                    '-Q', '33',
                    '-q', '20',
                    '-p', '90'
                ]
                subprocess.run(cmd, check=True)
                os.remove(intermediate_outname)
                logging.info(f"Quality filtered {intermediate_outname} to {final_out} and deleted intermediate file {intermediate_outname}")
            except Exception as e:
                logging.error(f"Failed to quality filter {intermediate_outname}: {e}")
                print(f"Error: Failed to quality filter {intermediate_outname}: {e}", file=sys.stderr)
                sys.exit(1)

    logging.info("Demultiplexing and quality filtering completed")


if __name__ == "__main__":
    main()



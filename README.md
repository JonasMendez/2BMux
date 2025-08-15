# 2BMux: 2b-RAD Demultiplexing and Sample Sheet Generation Toolkit

This repository contains a suite of scripts designed to process raw, multiplexed 2b-RAD sequencing data. The workflow allows you to demultiplex FASTQ files based on inline barcodes, generate the necessary sample mapping files, and prepare your data for downstream analysis tools like `ipyrad`.

## Overview

This toolkit includes three main components:

1.  **`2BMux.py`**: The core script that performs demultiplexing. It identifies reads based on a restriction site, extracts inline barcodes, assigns reads to the correct sample, performs quality filtering, and handles the specific orientation issues common in 2b-RAD libraries.
2.  **`PlateRelate.py`**: A helper script to create a properly formatted sample-to-barcode map (`-idmap` file). It takes a simple plate layout matrix and a list of sample names and generates the detailed CSV file required by the main demultiplexing script.
3.  **`concatenatepairs.sh`**: A utility shell script to concatenate paired-end (R1/R2) FASTQ files into a single file per sample. This is often required for downstream tools like `ipyrad` that use a single-end assembly workflow for 2b-RAD data.

## General Workflow

The recommended workflow is as follows:

1.  **Prepare Sample Information**: Create your input files:
    *   A CSV file listing your `SampleID` and a unique `Index` for each (e.g.`Plate1_Samples_input.csv`).
    *   A CSV file representing your physical plate layout, using the indices from the file above (e.g.`Plate1Matrix_input.csv`).

2.  **Generate the ID Map**: Use `PlateRelate.py` to convert your plate layout and sample list into the final `PlateX_IDs.csv` file. This file maps each sample to its specific plate, row, and barcode.

3.  **Demultiplex Raw Data**: Run the main `2BMux.py` script on your raw, multiplexed `R1` and `R2` FASTQ files, using the `PlateX_IDs.csv` file as the `-idmap` input. This will produce clean, demultiplexed FASTQ files for each sample.

4.  **(Optional) Merge Paired-End Files**: If your downstream analysis requires it, use the `concatenatepairs.sh` script to combine the demultiplexed `_R1_.fastq` and `_R2_.fastq` files into a single merged FASTQ file for each sample.

---

## 1. `PlateRelate.py` - Sample Sheet Generator

This script automates the creation of the sample-to-barcode mapping file required by the main demultiplexer. It translates a simple plate layout into the precise format needed for the 2BMux.py script -idmap input.

### Features

*   Parses a simple CSV matrix representing your 8-row plate.
*   Handles plates with a variable number of columns (up to 12).
*   Ignores empty cells or cells with non-numeric text (e.g., "Blank").
*   Maps the column position to a pre-defined 4-bp barcode.
*   Accepts command-line arguments for plate number and input files.
*   For multiple plates, run the script independently for each plate, using the proper plate number, and then concatenate output files for each and remove extra headers.
### Usage

```sh
python PlateRelate.py -plate <PLATE_NUM> -matrix <MATRIX_FILE> -ids <IDS_FILE>
```

**Arguments:**

*   `-plate`: (Required) An integer representing the plate number (e.g., `1`).
*   `-matrix`: (Required) Path to your plate layout CSV file. This file should have **no headers** and contain sample indices in an 8-row, up to 12-column format.
*   `-ids`: (Required) Path to your csv containing sample IDs and their corresponding index number. This file **must have headers** `SampleID` and `Index`.

### Example

**Input File 1: `Plate1_Samples_input.csv` (`-ids`)**
```csv
SampleID,Index
Sample_A,1
Sample_B,2
Sample_C,3
etc...
```

**Input File 2: `Plate1Matrix_input.csv` (`-matrix`)**
```csv
1,2,3
13,14,Blank,16
```

**Command:**
```sh
python PlateRelate.py -plate 1 -matrix plate1matrix.csv -ids samples.csv
```

**Output File: `Plate1_IDs.csv`**
```csv
Plate,Row,3illBC,SampleID
1,1,ACAC,Sample_A
1,1,GTCT,Sample_B
1,1,TGGT,Sample_C
```

---

## 2. `2BMux.py` - The Demultiplexer

This is the core script for processing raw 2b-RAD data. It reads paired-end FASTQ files, identifies the restriction site and inline barcode, assigns reads to samples, trims adaptors, and performs quality filtering.

### Features

*   **Flexible Filename Parsing**: Extracts plate, row, and read designator from filenames based on user-defined positions from underscore delimited filename elements.
*   **Regex-based Pattern Matching**: Uses regular expressions to identify the restriction site, barcode, and adaptor sequences, making the script adaptable to different library preparations.
*   **Orientation-Aware Demultiplexing**: Handles the forward and reverse orientations of the 2b-RAD fragment, which determines whether the inline barcode is the primary `3illBC` or its reverse complement `antiBC`.
*   **Specialized Modes**: Includes `--strict_orientation` to discard reads with conflicting site/barcode evidence and `--antiBC_only` for libraries where only the `antiBC` is reliable for demultiplexing.
*   **Quality Filtering**: Integrates with `fastq_quality_filter` (from the FASTX-Toolkit) to output quality filtered reads.
*   **Logging**: Generates a `demux_2bRAD.log` file to track progress and errors.
*   **Summary Reports**: Creates a `_barcode_summary.csv` for each input file, detailing read counts for each barcode and orientation.

### Dependencies

*   **FASTX-Toolkit**: This script requires the `fastq_quality_filter` command to be available in your system's PATH. You can install it via `apt-get install fastx-toolkit` (Debian/Ubuntu) or other package managers.

### Usage

The script should be run in the directory containing your raw, multiplexed `.fastq.gz` files.

```sh
python 2BMux.py -plate <P_POS> -row <R_POS> -read <D_POS> -idmap <ID_MAP_FILE> [OPTIONS]
```

**Required Arguments:**

*   `-plate <P_POS>`: The numeric position of the plate identifier in the underscore-separated filename (e.g., if `Plate1_...`, use `-plate 1`).
*   `-row <R_POS>`: The numeric position of the row identifier in the filename. (e.g., if `Plate1_Row1...`, use `-row 2`)
*   `-read <D_POS>`: The 1-based position of the read designator (`R1`/`R2`) in the filename. (e.g., if `Plate1_Row1_fastqinfo_R1.fastq.gz`, use `-read 4`)
*   `-idmap <ID_MAP_FILE>`: Path to the sample sheet created by `PlateRelate.py` (e.g., `Plate1_IDs.csv`).

**Optional Arguments:**
* In practice I have always used --antiBC_only and --allR1 for my specific workflows, but I leave the options here in case you find utility in excluding them.
*   `--strict_orientation`: If set, discards reads where the site orientation and barcode type are unexpected (e.g., a forward-oriented site with an `antiBC`). Default is to allow and flag them.
*   `--antiBC_only`: If set, only uses `antiBC` barcodes for demultiplexing, ignoring `3illBC` barcodes. This is specific for protocols using chimeric TCAC/GTGA barcodes which will require the use of a single barcode orientation for proper demultiplexing.
*   `--allR1`: If set, keeps demultiplexed R1 reads even if they do not have a corresponding R2 read.
*   `-min_bc_count <INT>`: Minimum number of reads required for a barcode to be written to an output file. Default: `10000`.
*   `-site <REGEX>`: The regex pattern for the restriction site. Defaults to the BcgI 2b-RAD pattern. i.e.: `'^(?:(?P<forward>.{12}CGA.{6}TGC.{12})|(?P<reverse>.{12}GCA.{6}TCG.{12}))'`
*   `-trim <INT>`: Number of bases to trim from both ends of the DNA fragment after demultiplexing. Default: `0`.

---

## 3. `concatenatepairs.sh` - FASTQ Merger

A simple shell script to concatenate R1 and R2 files into a single merged file, which is a common requirement for `ipyrad`'s `2brad` assembly mode.

### Features

*   Automatically finds all `_R1` files in a directory.
*   Handles both `.fastq` and `.fastq.gz` files.
*   Constructs the corresponding `_R2` and output filenames.
*   Does not delete original files.

### Usage

1.  Place the script in the directory containing your demultiplexed FASTQ files.
2.  Make it executable: `chmod +x concatenatepairs.sh`
3.  Run it: `bash concatenatepairs.sh`

### Example

**If your directory contains:**
*   `SampleA_R1_.fastq.gz`
*   `SampleA_R2_.fastq.gz`

**After running the script, you will have a new file:**
*   `SampleA.fastq.gz`

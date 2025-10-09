#!/usr/bin/env python3

import os
import re
import glob
import sqlite3
import csv
import argparse

# Function to extract the sample name from the contigs sequence name
def extract_sample_name_from_contig(contig_name):
    match = re.search(r'__(.*?)__', contig_name)
    return match.group(1) if match else "NA"

# Function to parse FASTA file and get header and sequences (handle multi-sequences per file)
def parse_fasta(fasta_path):
    sequences = []
    if not os.path.exists(fasta_path):
        return sequences
    with open(fasta_path) as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    sequences.append((header, ''.join(seq)))
                header = line
                seq = []
            else:
                seq.append(line)
        if header:
            sequences.append((header, ''.join(seq)))
    return sequences

# Function to extract the contig name from the header sequence
def extract_contig_name(header_line):
    match = re.search(r'contig:([^\|]+)', header_line)
    return match.group(1) if match else "NA"

# Function to get taxonomy from MAGNET.db
def get_contig_taxonomies(sqlite_db_path):
    taxo_dict = {}
    if not os.path.exists(sqlite_db_path):
        return taxo_dict
    try:
        conn = sqlite3.connect(sqlite_db_path)
        cursor = conn.cursor()
        cursor.execute("SELECT scaffold_name, gtdb_domain, gtdb_phylum, gtdb_class, gtdb_order, gtdb_family, gtdb_genus, gtdb_species FROM scaffolds;")
        for row in cursor.fetchall():
            scaffold_name = row[0]
            taxonomy = ";".join(row[1:])
            taxo_dict[scaffold_name] = taxonomy
    except Exception as e:
        print(f"[ERROR] SQLite read failed: {e}")
    finally:
        conn.close()
    return taxo_dict

# Function to extract bin_id frome header sequence
def extract_bin_id(header_line):
    match = re.search(r'bin_id:([^\|\s]+)', header_line)
    return match.group(1) if match else "NA"

# Function to get sequences bin by bin : CAUTION  with the path
def get_bin_sequences(sample_path):
    bin_seq_data = {}
    bin_files = glob.glob(os.path.join(sample_path, "MAGnet/ANVIO-SUMMARY/bin_by_bin/*/*-Ribosomal_RNA_16S-hmm-sequences.txt"))
    for file in bin_files:
        with open(file) as f:
            header, seq = None, []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header:
                        contig = extract_contig_name(header)
                        bin_id = extract_bin_id(header)
                        bin_seq_data[contig] = {
                            'bin_id': bin_id or "NA",
                            'header': header or "NA",
                            'sequence': ''.join(seq) or "NA"
                        }
                    header = line
                    seq = []
                else:
                    seq.append(line)
            # Last sequence
            if header:
                contig = extract_contig_name(header)
                bin_id = extract_bin_id(header)
                bin_seq_data[contig] = {
                    'bin_id': bin_id or "NA",
                    'header': header or "NA",
                    'sequence': ''.join(seq) or "NA"
                }
    return bin_seq_data

# Function to get bin taxonomy from GTDB-tk output (archaea and bacteria) : CAUTION with the path
def get_bin_taxonomies(sample_path):
    bin_taxo = {}
    summary_files = [
        os.path.join(sample_path, "MAGnet/GTDB-tk/output/gtdbtk.bac120.summary.tsv"),
        os.path.join(sample_path, "MAGnet/GTDB-tk/output/gtdbtk.ar53.summary.tsv")
    ]
    for file in summary_files:
        if os.path.exists(file):
            with open(file) as f:
                next(f)  # Skip header
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        bin_name, taxonomy = parts[0], parts[1]
                        bin_taxo[bin_name] = taxonomy
    return bin_taxo

# Function to ensure space and columns informations : replace with NA
def safe(value):
    return value if value else "NA"

# Function to process each sample from a list and create .tsv file : CAUTION with path
def process_sample(sample_name):
    rows = []
    sample_path = sample_name
    fasta_path = os.path.join(sample_path, "assembly/datatables/Ribosomal_RNA_16S.fa")
    sqlite_path = os.path.join(sample_path, "MAGnet/MAGNET.db")

    contig_taxonomies = get_contig_taxonomies(sqlite_path)
    bin_sequences = get_bin_sequences(sample_path)
    bin_taxonomies = get_bin_taxonomies(sample_path)

    sequences = parse_fasta(fasta_path)
    for header, sequence in sequences:
        contig = extract_contig_name(header)
        sample_from_contig = extract_sample_name_from_contig(contig)
        contig_taxo = contig_taxonomies.get(contig, "NA")

        bin_info = bin_sequences.get(contig, {})
        bin_id = safe(bin_info.get("bin_id"))
        bin_taxo = bin_taxonomies.get(bin_id, "NA")
        bin_header = safe(bin_info.get("header"))
        bin_sequence = safe(bin_info.get("sequence"))

        # Duplicate sequence name in sequence columns, formatted as FASTA
        seq_assembly = f"{header}\n{sequence}"
        seq_bin = f"{bin_header}\n{bin_sequence}" if bin_sequence != "NA" else "NA"

        row = list(map(safe, [
            sample_from_contig,    # Column 1: sample name
            header,                # Column 2: sequence header assembly
            seq_assembly,          # Column 3: sequence assembly with header
            contig_taxo,           # Column 4: contig taxonomy
            bin_id,                # Column 5: bin id
            bin_taxo,              # Column 6: bin taxonomy
            bin_header,            # Column 7: sequence header bin
            seq_bin                # Column 8: sequence bin with header
        ]))
        rows.append(row)
    return rows

def main():
    parser = argparse.ArgumentParser(description="Extract 16S sequence info from samples")
    parser.add_argument("--input", "-i", required=True, help="Path to file listing samples (one per line)")
    parser.add_argument("--output", "-o", required=True, help="Path to output TSV file")
    args = parser.parse_args()

    with open(args.input) as f:
        samples = [line.strip() for line in f if line.strip()]

    print(f"[INFO] Samples to process: {samples}")

    all_rows = []
    for sample in samples:
        print(f"[INFO] Processing sample: {sample}")
        rows = process_sample(sample)
        print(f"[INFO] Found {len(rows)} sequences in sample {sample}")
        all_rows.extend(rows)

    with open(args.output, "w", newline='') as f:
        writer = csv.writer(f, delimiter="\t", quoting=csv.QUOTE_MINIMAL)
        writer.writerow([
            "Sample",
            "SeqName_Assembly",
            "Sequence_Assembly",
            "Contig_Taxonomy",
            "Bin_ID",
            "Bin_Taxonomy",
            "SeqName_Bin",
            "Sequence_Bin"
        ])
        writer.writerows(all_rows)

    print(f"[DONE] Output written to: {args.output}")

if __name__ == "__main__":
    main()

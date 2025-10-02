#!/usr/bin/env python3

import os
import gzip
import argparse

# Function to handle .faa and .faa.gz file
def open_faa(filepath):
    return gzip.open(filepath, 'rt') if filepath.endswith('.gz') else open(filepath, 'r')

# Function to parse sequences headers and return parsed info
def parse_header(header_line, genome_prefix='NA'):
    header_line = header_line.strip().lstrip('>')
    if '|' in header_line:
        _, rest = header_line.split('|', 1)  # ignore genome name
    else:
        rest = header_line

    parts = rest.split('#')
    if len(parts) < 4:
        return None
    full_id = parts[0].strip()
    start = parts[1].strip()
    stop = parts[2].strip()
    strand = parts[3].strip()

    if '_' in full_id:
        scaffold = '_'.join(full_id.split('_')[:-1])
        orf_number = full_id.split('_')[-1]
    else:
        scaffold = full_id
        orf_number = 'NA'

    unique_orf_id = f"{scaffold}_{orf_number}"
    return [unique_orf_id, genome_prefix, scaffold, start, stop, strand]

# Parse headers for a specific .faa file and genome name
def parse_faa_file(filepath, genome_name):
    results = []
    with open_faa(filepath) as f:
        for line in f:
            if line.startswith('>'):
                parsed = parse_header(line, genome_prefix=genome_name)
                if parsed:
                    results.append(parsed)
    return results

# Search protein.faa files
def collect_faa_files(base_dir):
    faa_files = []
    for root, _, files in os.walk(base_dir):
        for file in files:
            if file.endswith('protein.faa') or file.endswith('protein.faa.gz') or file.endswith('proteins.faa'):
                faa_files.append(os.path.join(root, file))
    return faa_files

# Concat .faa and rename headers with unique ORF ID only
def concat_faa_files(faa_files, output_fasta):
    with open(output_fasta, 'w') as out_f:
        for filepath in faa_files:
            with open_faa(filepath) as in_f:
                for line in in_f:
                    if line.startswith('>'):
                        parsed = parse_header(line)
                        if parsed:
                            unique_orf_id = parsed[0]
                            out_f.write(f">{unique_orf_id}\n")
                    else:
                        out_f.write(line)
    print(f"All files concatenated into: {output_fasta}")

def main():
    parser = argparse.ArgumentParser(description="Concatenate and parse protein.faa files from GTDB-Tk marker genes")
    parser.add_argument("-i", "--input", help="Base directory containing genome subfolders with protein.faa files")
    parser.add_argument("-o", "--output", default="dataset_feature.tsv", help="Output TSV filename for parsed headers")
    parser.add_argument("-c", "--concat_output", default="dataset_protein.faa", help="Output concatenated fasta filename")
    parser.add_argument("-g", "--genome_name_from_file", action='store_true', help="Use the filename (without suffix) as genome name instead of parent folder")
    args = parser.parse_args()

    faa_files = collect_faa_files(args.input)
    if not faa_files:
        print("No protein.faa files found in the specified directory.")
        return

    concat_faa_files(faa_files, args.concat_output)

    results = []

    for filepath in faa_files:
        if args.genome_name_from_file:
            filename = os.path.basename(filepath)
            if filename.endswith('_proteins.faa'):
                genome_name = filename[:-len('_proteins.faa')]
            elif filename.endswith('_proteins.faa.gz'):
                genome_name = filename[:-len('_proteins.faa.gz')]
            elif filename.endswith('_protein.faa'):
                genome_name = filename[:-len('_protein.faa')]
            else:
                genome_name = os.path.splitext(filename)[0]
        else:
            genome_name = os.path.basename(os.path.dirname(filepath))
        results.extend(parse_faa_file(filepath, genome_name))

    with open(args.output, 'w') as out:
        out.write("orf\tgenome\tscaffold\tstart\tstop\tstrand\n")
        for row in results:
            out.write('\t'.join(row) + '\n')

    print(f"Parsed TSV file created: {args.output} ({len(results)} ORFs)")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import random
import gzip
import shutil
import re
from urllib.request import urlretrieve
from pathlib import Path

GTDB_BASE_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/"
ARCHAEA_META = "ar53_metadata.tsv.gz"
BACTERIA_META = "bac120_metadata.tsv.gz"

# Function to clean taxon name
def sanitize_taxon_name(name):
    name = name.strip()
    name = re.sub(r'[^\w\-]', '_', name)
    name = re.sub(r'_+', '_', name)
    return name

# Function to choose domain_file, create meta_files
def download_metadata(domain, outdir):
    filename_gz = ARCHAEA_META if domain == "Archaea" else BACTERIA_META
    filename_tsv = filename_gz[:-3]
    parent_dir = os.path.dirname(os.path.abspath(outdir))
    out_path_tsv = os.path.join(parent_dir, filename_tsv)
    out_path_gz = os.path.join(parent_dir, filename_gz)

    if os.path.exists(out_path_tsv):
        print(f"Metadata file already exists: {out_path_tsv}")
        return out_path_tsv

    if not os.path.exists(out_path_gz):
        url = os.path.join(GTDB_BASE_URL, filename_gz)
        print(f"Downloading metadata archive: {url}")
        urlretrieve(url, out_path_gz)
    else:
        print(f"Metadata archive already downloaded: {out_path_gz}")

    print(f"Decompressing {out_path_gz} to {out_path_tsv}")
    with gzip.open(out_path_gz, 'rb') as f_in, open(out_path_tsv, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

    return out_path_tsv

# Function to load meta_file with pandas
def load_metadata(file_path):
    return pd.read_csv(file_path, sep='\t', low_memory=False)

# Function to filter genome by taxonomy levels ("g","Desulfamplus" = g__Desulfamplus)
def filter_by_taxonomy(df, tax_level, tax_name):
    pattern = rf"\b{tax_level}__{tax_name}\b"
    return df[df['gtdb_taxonomy'].str.contains(pattern, regex=True, na=False)]

# Function to keep only representative genomes
def filter_representative(df):
    return df[df['gtdb_representative'] == 't']

# Function to extract the last gtdb taxon level (to rename assembly), keeping rank prefix
def extract_most_specific_taxon(taxonomy):
    parts = taxonomy.split(';')
    for part in reversed(parts):
        if '__' in part and len(part) > 3:
            return sanitize_taxon_name(part)  # ex: s__coli -> s_coli
    return "unknown_taxon"

# Function to create link and download genomes
def download_and_process_genome(accession, assembly_name, taxonomy, outdir):
    prefix = accession.split('_')[0]
    accession_numbers = accession.split('_')[1].split('.')[0]
    part1 = accession_numbers[0:3]
    part2 = accession_numbers[3:6]
    part3 = accession_numbers[6:9]
    safe_assembly_name = assembly_name.replace(' ', '_')
    folder_name = f"{accession}_{safe_assembly_name}"
    
    fasta_url = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{part1}/{part2}/{part3}/{folder_name}/{folder_name}_genomic.fna.gz"
    temp_gz_path = os.path.join(outdir, f"{accession}.temp.fna.gz")

    # Final filename based on taxonomy
    final_taxon_name = extract_most_specific_taxon(taxonomy)
    filename = f"{final_taxon_name}.fna"
    final_path = os.path.join(outdir, filename)

    # Ensure unique filename
    counter = 1
    while os.path.exists(final_path):
        filename = f"{final_taxon_name}_{counter}.fna"
        final_path = os.path.join(outdir, filename)
        counter += 1

    try:
        print(f"Downloading: {fasta_url}")
        urlretrieve(fasta_url, temp_gz_path)

        # Rename headers in fasta: keep old assembly name in contig headers
        with gzip.open(temp_gz_path, 'rt') as f_in, open(final_path, 'w') as f_out:
            for line in f_in:
                if line.startswith('>'):
                    contig_name = line.strip()[1:]
                    new_header = f">{accession}|{contig_name}\n"
                    f_out.write(new_header)
                else:
                    f_out.write(line)
        os.remove(temp_gz_path)
        return filename
    except Exception as e:
        print(f"Download or processing failed for {accession}: {e}")
        return None

# Function to define some filters and process genomes
def process_filters(df, filters, max_per_filter, output_dir, representative_only, global_tsv_list, code_label):
    translation_table = code_label.strip('_')
    for tax_filter in filters:
        if '__' not in tax_filter:
            print(f"Invalid filter: {tax_filter}")
            continue
        level, name = tax_filter.split('__', 1)
        filtered_df = filter_by_taxonomy(df, level, name)
        if representative_only:
            filtered_df = filter_representative(filtered_df)
        if filtered_df.empty:
            print(f"No genomes found for filter {tax_filter}")
            continue

        # Extract required columns and random sampling
        filtered_accessions = filtered_df[['ncbi_genbank_assembly_accession', 'ncbi_assembly_name', 'gtdb_taxonomy']].dropna(subset=['ncbi_genbank_assembly_accession', 'ncbi_assembly_name', 'gtdb_taxonomy'])
        if max_per_filter is not None and len(filtered_accessions) > max_per_filter:
            filtered_accessions = filtered_accessions.sample(n=max_per_filter, random_state=42)

        print(f"{len(filtered_accessions)} genomes selected for filter {tax_filter}.")

        # Download and process each genome
        for _, row in filtered_accessions.iterrows():
            acc = row['ncbi_genbank_assembly_accession']
            asm_name = row['ncbi_assembly_name']
            taxo = row['gtdb_taxonomy']

            filename = download_and_process_genome(acc, asm_name, taxo, output_dir)
            if filename:
                global_tsv_list.append((acc, taxo, translation_table, filename))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--domain", required=True, choices=["Bacteria", "Archaea"])
    parser.add_argument("--tax_filter", action="append", default=[], help="Main group filters to store in code11")
    parser.add_argument("--tax_filter_code25", action="append", default=[], help="Main group filters to store in code25")
    parser.add_argument("--outgroup_filter_code11", action="append", default=[], help="Outgroup filters to store in code11 (representative only)")
    parser.add_argument("--outgroup_filter_code25", action="append", default=[], help="Outgroup filters to store in code25 (representative only)")
    parser.add_argument("--max_main_per_filter", type=int, default=None)
    parser.add_argument("--max_main_per_filter_code25", type=int, default=None)
    parser.add_argument("--max_outgroup_per_filter", type=int, default=None)
    parser.add_argument("--outdir", default="genomes", help="Base output directory")

    args = parser.parse_args()

    random.seed(42)

    Path(args.outdir).mkdir(parents=True, exist_ok=True)
    genomes_dir = os.path.join(args.outdir, "genomes")

    # Download and load metadata
    meta_file = download_metadata(args.domain, args.outdir)
    df = load_metadata(meta_file)
    all_tsv_entries = []

    if args.tax_filter:
        code11_dir = os.path.join(genomes_dir, "code11")
        Path(code11_dir).mkdir(parents=True, exist_ok=True)
        process_filters(df, args.tax_filter, args.max_main_per_filter, code11_dir, False, all_tsv_entries, "_code11")

    if args.tax_filter_code25:
        code25_dir = os.path.join(genomes_dir, "code25")
        Path(code25_dir).mkdir(parents=True, exist_ok=True)
        process_filters(df, args.tax_filter_code25, args.max_main_per_filter_code25, code25_dir, False, all_tsv_entries, "_code25")

    if args.outgroup_filter_code11:
        outgroup_dir11 = os.path.join(genomes_dir, "code11")
        Path(outgroup_dir11).mkdir(parents=True, exist_ok=True)
        process_filters(df, args.outgroup_filter_code11, args.max_outgroup_per_filter, outgroup_dir11, True, all_tsv_entries, "_outgroup_code11")

    if args.outgroup_filter_code25:
        outgroup_dir25 = os.path.join(genomes_dir, "code25")
        Path(outgroup_dir25).mkdir(parents=True, exist_ok=True)
        process_filters(df, args.outgroup_filter_code25, args.max_outgroup_per_filter, outgroup_dir25, True, all_tsv_entries, "_outgroup_code25")

    # Create metadata summary file if entries exist
    if all_tsv_entries:
        final_tsv_path = os.path.join(args.outdir, "genomes.metadata.txt")
        metadata_df = pd.DataFrame(all_tsv_entries, columns=["assembly_accession", "taxonomy", "translation_table", "filename"])
        metadata_df.to_csv(final_tsv_path, sep='\t', index=False)
        print(f"Metadata file written to: {final_tsv_path}")

if __name__ == "__main__":
    main()

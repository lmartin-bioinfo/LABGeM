#!/usr/bin/env python3

import argparse
import pandas as pd
import os

# Function to extract taxonomy by taxonomic level
def extract_taxon_at_level(taxonomy: str, level: str) -> str:
    if pd.isna(taxonomy):
        return "NA"
    parts = taxonomy.split(';')
    prefix = level + "__"
    for part in parts:
        if part.startswith(prefix):
            return part[len(prefix):]
    return "NA"

def main():
    parser = argparse.ArgumentParser(description="Create genome-to-taxonomy files for iTOL")
    parser.add_argument("--metadata_file", required=True, help="Path to genomes.metadata.txt")
    parser.add_argument("--summary_file", required=True, help="Path to gtdbtk.bac120.summary.tsv")
    parser.add_argument("--levels", required=True, nargs='+', choices=['d', 'p', 'c', 'o', 'f', 'g', 's'], help="Taxonomic levels to extract")
    parser.add_argument("--output_dir", default=".", help="Output directory for result files")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load metadata
    metadata = pd.read_csv(args.metadata_file, sep='\t', dtype=str)
    metadata['user_genome'] = metadata['filename'].str.replace('.fna$', '', regex=True)

    # Load GTDB-Tk output
    gtdbtk = pd.read_csv(args.summary_file, sep='\t', dtype=str)
    gtdbtk = gtdbtk[['user_genome', 'classification']].rename(columns={'classification': 'taxonomy'})

    # Subset metadata
    metadata = metadata[['user_genome', 'taxonomy']]

    # Identify new genomes (not present in metadata)
    existing_user_genomes = set(metadata['user_genome'])
    gtdbtk_new = gtdbtk[~gtdbtk['user_genome'].isin(existing_user_genomes)]

    # Combine data: priority to genomes.metadata.txt taxonomy (latest taxonomy from gtdb)
    combined = pd.concat([metadata, gtdbtk_new], ignore_index=True)

    # Generate taxonomy files for each level
    for level in args.levels:
        combined[f'taxonomy_{level}'] = combined['taxonomy'].apply(lambda x: extract_taxon_at_level(x, level))
        output_file = os.path.join(args.output_dir, f'genome2taxonomy__{level}.txt')
        output_df = combined[['user_genome', f'taxonomy_{level}']].rename(
            columns={'user_genome': 'assembly', f'taxonomy_{level}': 'taxonomy'}
        )
        output_df.to_csv(output_file, sep='\t', index=False)
        print(f"Saved: {output_file}")

if __name__ == "__main__":
    main()

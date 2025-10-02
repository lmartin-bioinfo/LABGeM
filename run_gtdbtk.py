#!/usr/bin/env python

import argparse
import gzip
import shutil
import subprocess
import re
import os
from pathlib import Path
from Bio import SeqIO

# Function to check if a FASTA file is valid by parsing sequences
def is_valid_fasta(fasta_path: Path) -> bool:
    try:
        with open(fasta_path, "r") as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            return len(records) > 0
    except Exception as e:
        print(f" Error in {fasta_path.name}: {e}")
        return False

# Function to decompress .gz FASTA files to plain text
def decompress_fasta_gz(src_path: Path, dst_path: Path):
    with gzip.open(src_path, 'rt') as f_in, open(dst_path, 'w') as f_out:
        shutil.copyfileobj(f_in, f_out)

# Function to process each genome file: decompress if needed and rename to .fna
def process_genome_file(file_path: Path, output_dir: Path) -> bool:
    # Extract base name (without extension and .gz if present)
    if file_path.suffix == ".gz":
        base_name = Path(file_path.stem).stem
    else:
        base_name = file_path.stem

    dst_path = output_dir / f"{base_name}.fna"

    if dst_path.exists():
        print(f"{dst_path.name} already exists. Skipping.")
        return False

    if file_path.suffix == ".gz":
        temp_path = output_dir / f".temp_{base_name}.fna"
        decompress_fasta_gz(file_path, temp_path)
        if is_valid_fasta(temp_path):
            temp_path.rename(dst_path)
            print(f"Valid and decompressed: {dst_path.name}")
            return True
        else:
            print(f"Invalid FASTA content (gz): {file_path.name}")
            temp_path.unlink()
            return False

    elif file_path.suffix in [".fasta", ".fna"]:
        if is_valid_fasta(file_path):
            shutil.copy2(file_path, dst_path)
            print(f"Valid FASTA copied and renamed to: {dst_path.name}")
            return True
        else:
            print(f"Invalid FASTA file: {file_path.name}")
            return False

    else:
        print(f"Skipped (not .fasta, .fna, or .gz): {file_path.name}")
        return False

# Function to walk through input directories and process genome files
def collect_genomes(input_dirs, output_dir: Path, filtered_dir: list = None, name_filter: str = ""):
    output_dir.mkdir(parents=True, exist_ok=True)
    filtered_dirs_resolved = set(Path(fd).resolve() for fd in filtered_dir) if filtered_dir else set()

    for dir_path in input_dirs:
        dir_path = Path(dir_path)
        print(f"\nProcessing directory: {dir_path}")
        for file in dir_path.rglob("*"):
            if not file.is_file():
                continue

            # Apply name filtering only on specified directories if needed
            is_filtered = dir_path.resolve() in filtered_dirs_resolved
            if is_filtered and name_filter and not any(keyword in file.name for keyword in name_filter):
                continue

            process_genome_file(file, output_dir)

# Function to run GTDB-Tk classify workflow
def run_gtdbtk_classify(genome_dir: Path, output_dir: Path, cpus: int = 4, extension: str = "fasta"):
    print("\nRunning GTDB-Tk classify_wf with --genome_dir...\n")
    cmd = [
        "gtdbtk", "classify_wf",
        "--genome_dir", str(genome_dir),
        "--out_dir", str(output_dir),
        "--cpus", str(cpus),
        "--extension", extension,
        "--skip_ani_screen",
        "--keep_intermediates"
    ]
    subprocess.run(cmd, check=True)
    print(f"\nGTDB-Tk classify_wf finished. Results at: {output_dir}")


# Function to run GTDB-Tk classify workflow with a batchfile (you can specify genetic code)
def run_gtdbtk_classify_batch(batchfile: Path, output_dir: Path, cpus: int = 4, extension: str = "fasta"):
    print("\nRunning GTDB-Tk classify_wf with --batchfile...\n")
    cmd = [
        "gtdbtk", "classify_wf",
        "--batchfile", str(batchfile),
        "--out_dir", str(output_dir),
        "--cpus", str(cpus),
        "--extension", extension,
        "--skip_ani_screen",
        "--keep_intermediates"
    ]
    subprocess.run(cmd, check=True)
    print(f"\nGTDB-Tk classify_wf finished. Results at: {output_dir}")

# Function to run GTDB-Tk de_novo_wf (phylogenetic tree construction)
def run_gtdbtk_de_novo(genome_dir: Path, output_dir: Path, outgroup_taxon: str, domain: str = "bacteria", cpus: int = 4, extension: str = "fasta"):
    if domain not in ["bacteria", "archaea"]:
        raise ValueError("Domain must be 'bacteria' or 'archaea'")
    print("\nRunning GTDB-Tk de_novo workflow for tree construction...\n")
    cmd = [
        "gtdbtk", "de_novo_wf",
        "--genome_dir", str(genome_dir),
        f"--{domain}",
        "--outgroup_taxon", outgroup_taxon,
        "--out_dir", str(output_dir),
        "--cpus", str(cpus),
        "--extension", extension
    ]
    subprocess.run(cmd, check=True)
    print(f"\nGTDB-Tk de_novo workflow finished. Tree results at: {output_dir}")

def main():
    parser = argparse.ArgumentParser(description="Prepare and run GTDB-Tk on genome files.")
    parser.add_argument("--genome_dirs", nargs="+", help="List of directories containing genome files")
    parser.add_argument("--filtered_dir", nargs="+", help="Directories to apply filename filtering to")
    parser.add_argument("--filter_name", nargs="+", default="", help="Filter filenames by keyword (only in filtered_dir)")
    parser.add_argument("--tmp_input_dir", default="input_gtdbtk", help="Temporary input directory for GTDB-Tk")
    parser.add_argument("--output_dir", default="output_gtdbtk", help="Output directory for GTDB-Tk results")
    parser.add_argument("--cpus", type=int, default=4, help="Number of CPUs for GTDB-Tk")
    parser.add_argument("--extension", default="fasta", help="Genome file extension (default: 'fasta')")
    parser.add_argument("--run_classify_wf", action="store_true", help="Run GTDB-Tk classify workflow")
    parser.add_argument("--de_novo_tree", action="store_true", help="Run GTDB-Tk de_novo workflow to build phylogenetic tree")
    parser.add_argument("--outgroup_taxon", type=str, default=None, help="GTDB taxon name to use as outgroup (e.g., 'p__Firmicutes')")
    parser.add_argument("--domain", type=str, default="bacteria", choices=["bacteria", "archaea"], help="Domain of the genomes (bacteria or archaea)")
    parser.add_argument("--batchfile", type=str, help="Optional batchfile for GTDB-Tk classify_wf. Tab-separated file with path/to/genome \t genome_name \t translation_table")

    args = parser.parse_args()

    output = Path(args.output_dir)

    # If no batchfile, process genome_dirs
    if not args.batchfile:
        if not args.genome_dirs:
            raise ValueError("You must provide either --genome_dirs or --batchfile.")
        tmp_input = Path(args.tmp_input_dir)
        tmp_input.mkdir(parents=True, exist_ok=True)

        collect_genomes(
            input_dirs=args.genome_dirs,
            output_dir=tmp_input,
            filtered_dir=args.filtered_dir,
            name_filter=args.filter_name
        )
    else:
        tmp_input = None  # not used with batchfile

    # Run classify_wf if requested
    if args.run_classify_wf:
        if args.batchfile:
            run_gtdbtk_classify_batch(batchfile=Path(args.batchfile), output_dir=output, cpus=args.cpus, extension=args.extension)
        else:
            run_gtdbtk_classify(genome_dir=tmp_input, output_dir=output, cpus=args.cpus, extension=args.extension)

    # Run de_novo_wf if requested
    if args.de_novo_tree:
        if not args.outgroup_taxon:
            raise ValueError("You must specify --outgroup_taxon to run de_novo_wf.")
        run_gtdbtk_de_novo(genome_dir=tmp_input, output_dir=output, outgroup_taxon=args.outgroup_taxon, domain=args.domain, cpus=args.cpus, extension=args.extension)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
animalfinders.py — Find the top organism match for a DNA/protein sequence via BLAST,
and return its Wikipedia page link.

Usage:
    Just run:
        ./animalfinders.py

Requirements:
    - Biopython (install via `pip install biopython`)
    - Internet connection (BLAST query to NCBI)
"""

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import urllib.parse
import sys


def run_blast(sequence: str, program: str = "blastn", database: str = "nt"):
    """Run BLAST search on NCBI and return parsed XML results."""
    print("[INFO] Running BLAST search... this may take a moment.")
    result_handle = NCBIWWW.qblast(program, database, sequence)
    print("[INFO] BLAST search completed. Parsing results...")
    blast_record = NCBIXML.read(result_handle)
    return blast_record


def extract_top_hit(blast_record):
    """Return organism name from the top BLAST hit."""
    if not blast_record.alignments:
        print("[WARN] No hits found.")
        return None
    title = blast_record.alignments[0].title
    # Example title: "gi|12345|ref|NM_12345.1| Homo sapiens gene ABC"
    # Extract organism name heuristically
    parts = title.split()
    for i in range(len(parts) - 1):
        if parts[i][0].isupper() and parts[i+1][0].islower():
            organism = f"{parts[i]} {parts[i+1]}"
            return organism
    return None


def main():
    print("=== Animal Finder ===")
    print("1) Enter sequence manually")
    print("2) Provide FASTA file path")
    choice = input("Choose input method (1 or 2): ").strip()

    if choice == "1":
        sequence = input("Enter your DNA/protein sequence: ").strip()
    elif choice == "2":
        fasta_path = input("Enter path to FASTA file: ").strip()
        try:
            record = SeqIO.read(fasta_path, "fasta")
            sequence = str(record.seq)
        except Exception as e:
            sys.exit(f"[ERROR] Failed to read FASTA file: {e}")
    else:
        sys.exit("[ERROR] Invalid choice. Please select 1 or 2.")

    # Run BLAST
    blast_record = run_blast(sequence)
    organism = extract_top_hit(blast_record)

    if not organism:
        print("[INFO] Could not identify organism from top hit.")
        sys.exit(0)

    # Create Wikipedia link
    wiki_name = urllib.parse.quote(organism)
    wiki_link = f"https://en.wikipedia.org/wiki/{wiki_name}"

    print(f"\n[RESULT] Top organism match: {organism}")
    print(f"[LINK] {wiki_link}\n")


if __name__ == "__main__":
    main()


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

def calculate_gc(sequence: str) -> float:
    """Return GC content percentage of a sequence."""
    g = sequence.upper().count("G")
    c = sequence.upper().count("C")
    return round((g + c) / len(sequence) * 100, 2) if len(sequence) > 0 else 0.0

def get_sequence_type(sequence: str) -> str:
    """Guess whether the sequence is DNA or protein."""
    dna_letters = set("ACGTU")
    if set(sequence.upper()) <= dna_letters:
        return "DNA"
    else:
        return "Protein"

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

    # Determine sequence type and length
    seq_type = get_sequence_type(sequence)
    seq_length = len(sequence)
    gc_content = calculate_gc(sequence)

    # Run BLAST
    blast_record = run_blast(sequence)

    if not blast_record.alignments:
        print("[INFO] No BLAST hits found.")
        sys.exit(0)

    # Get top hit full scientific name
    top_hit = blast_record.alignments[0].hit_def
    wiki_name = urllib.parse.quote(top_hit)
    wiki_link = f"https://en.wikipedia.org/wiki/{wiki_name}"

    # Display results
    print("\n=== Results ===")
    print(f"Full scientific name (top hit): {top_hit}")
    print(f"Wikipedia link: {wiki_link}")
    print(f"Sequence type: {seq_type}")
    print(f"Sequence length: {seq_length} bp/amino acids")
    print(f"GC content: {gc_content}%\n")

if __name__ == "__main__":
    main()



#!/usr/bin/env python3
"""
animalfinders.py — Find the top organism match for a DNA/protein sequence via BLAST,
and return its Wikipedia page link.

Usage:
    ./animalfinders.py

Requirements:
    - Biopython (`pip install biopython`)
    - Requests (`pip install requests`)
    - Internet connection (BLAST query to NCBI)
"""

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import urllib.parse
import requests
import sys


# --------------------------
# Utility Functions
# --------------------------

def run_blast(sequence: str, program: str = "blastn", database: str = "nt"):
    """Run BLAST search on NCBI and return parsed XML results."""
    print(f"[INFO] Running {program.upper()} search... this may take a while.")
    try:
        result_handle = NCBIWWW.qblast(program, database, sequence)
        print("[INFO] BLAST search completed. Parsing results...")
        blast_record = NCBIXML.read(result_handle)
        return blast_record
    except Exception as e:
        sys.exit(f"[ERROR] BLAST query failed: {e}")


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
    return "Protein"


def clean_hit_name(hit_def: str) -> str:
    """Extract the most likely scientific name from a BLAST hit description."""
    # Remove extra biological terms
    keywords = [
        "mitochondrion", "chloroplast", "plasmid", "chromosome",
        "complete genome", "genome assembly", "isolate", "strain",
        "DNA", "RNA", "sequence"
    ]
    for word in keywords:
        if word.lower() in hit_def.lower():
            hit_def = hit_def.split(word, 1)[0]
            break

    # Clean up punctuation and parentheses
    clean_name = hit_def.replace(",", "").split("(")[0].strip()

    # Keep only the genus + species if possible
    parts = clean_name.split()
    if len(parts) >= 2:
        clean_name = " ".join(parts[:2])

    return clean_name


def get_wiki_link(scientific_name: str) -> str:
    """Return a working Wikipedia link if available, otherwise a Google search."""
    base_url = "https://en.wikipedia.org/wiki/"
    encoded = urllib.parse.quote(scientific_name.replace(" ", "_"))
    wiki_url = base_url + encoded

    try:
        response = requests.get(wiki_url, timeout=5)
        if response.status_code == 200:
            return wiki_url
        else:
            raise ValueError("Not found")
    except Exception:
        query = urllib.parse.quote(scientific_name)
        return f"https://www.google.com/search?q={query}+site:en.wikipedia.org"


# --------------------------
# Main Function
# --------------------------

def main():
    print("=== 🧬 Animal Finder ===")
    print("1) Enter sequence manually")
    print("2) Provide FASTA file path")
    choice = input("Choose input method (1 or 2): ").strip()

    # Input handling
    if choice == "1":
        sequence = input("Enter your DN


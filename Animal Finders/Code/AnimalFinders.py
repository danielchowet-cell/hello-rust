#!/usr/bin/env python3
"""
Requirements:
    - Internet connection (for BLAST query)
    - Python 3
    - Biopython, Requests (auto-installed if missing)
"""

import sys
import subprocess
import importlib


# Auto-install dependencies


def ensure_package(pkg_name):
    """Ensure a package is installed, install it if not."""
    try:
        importlib.import_module(pkg_name)
    except ImportError:
        print(f"[INFO] Installing missing dependency: {pkg_name}")
        subprocess.check_call([sys.executable, "-m", "pip", "install", pkg_name])

for package in ["biopython", "requests"]:
    ensure_package(package)

# Importing
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import urllib.parse
import requests



# Utility Functions


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
    return "DNA" if set(sequence.upper()) <= dna_letters else "Protein"


def clean_hit_name(hit_def: str) -> str:
    """Extract the most likely scientific name from a BLAST hit description."""
    keywords = [
        "mitochondrion", "chloroplast", "plasmid", "chromosome",
        "complete genome", "genome assembly", "isolate", "strain",
        "DNA", "RNA", "sequence"
    ]
    for word in keywords:
        if word.lower() in hit_def.lower():
            hit_def = hit_def.split(word, 1)[0]
            break

    clean_name = hit_def.replace(",", "").split("(")[0].strip()
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


def save_results(filename: str, info: dict):
    """Save the BLAST summary results to a text file."""
    try:
        with open(filename, "w", encoding="utf-8") as f:
            f.write("=== Animal Finder Results ===\n\n")
            for key, value in info.items():
                f.write(f"{key}: {value}\n")
        print(f"\n[INFO] Results saved to '{filename}'")
    except Exception as e:
        print(f"[ERROR] Could not save results to file: {e}")



# Main Function


def main():
    print("=== 🧬 Animal Finder ===")
    print("1) Enter sequence manually")
    print("2) Provide FASTA file path")
    choice = input("Choose input method (1 or 2): ").strip()

    # Input handling
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

    # Basic sequence info
    seq_type = get_sequence_type(sequence)
    seq_length = len(sequence)
    gc_content = calculate_gc(sequence) if seq_type == "DNA" else "N/A"

    print(f"\n[INFO] Sequence type detected: {seq_type}")
    print(f"[INFO] Sequence length: {seq_length}")
    if seq_type == "DNA":
        print(f"[INFO] GC content: {gc_content}%")

    # Run BLAST
    program = "blastn" if seq_type == "DNA" else "blastp"
    blast_record = run_blast(sequence, program=program)

    # Check for results
    if not blast_record.alignments:
        print("[INFO] No BLAST hits found.")
        sys.exit(0)

    # Extract and clean hit name
    top_hit = blast_record.alignments[0].hit_def
    clean_name = clean_hit_name(top_hit)
    wiki_link = get_wiki_link(clean_name)

    # Prepare results dictionary
    results = {
        "Top BLAST hit": top_hit,
        "Most likely organism": clean_name,
        "Wikipedia link": wiki_link,
        "Sequence type": seq_type,
        "Sequence length": f"{seq_length} bp/amino acids",
        "GC content": f"{gc_content}%" if seq_type == "DNA" else "N/A"
    }

    # Display results
    print("\n=== 🔎 Results ===")
    for key, value in results.items():
        print(f"{key}: {value}")

    # Save to file
    save_results("results.txt", results)


if __name__ == "__main__":
    main()


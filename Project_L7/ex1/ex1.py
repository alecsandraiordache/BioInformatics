# -*- coding: utf-8 -*-
"""
Detects repetitions between 3 and 6 bases in a biological DNA sequence
(downloaded from NCBI, length 1000–3000 nt).
A repetition = a motif that occurs N or more times (overlaps allowed).
Results are printed and saved to a text file.
"""

from collections import defaultdict
import re

# === SETTINGS ===============================================================
FASTA_FILE = "escherichia_coli.fasta"   # NCBI DNA sequence (e.g. 16S rRNA)
K_MIN, K_MAX = 3, 6                     # motif lengths
MIN_REPEATS = 2                         # N times (minimum repeat count)
OUT_FILE = "repeats_3_6.txt"            # output text file
# ============================================================================

def read_fasta(path):
    """Read a FASTA file and return (header, cleaned DNA sequence)."""
    with open(path, "r") as f:
        lines = [l.strip() for l in f if l.strip()]
    header = lines[0] if lines[0].startswith(">") else ">unknown"
    seq = "".join(l for l in lines if not l.startswith(">")).upper()
    seq = seq.replace("U", "T")
    seq = re.sub(r"[^ACGT]", "", seq)
    return header, seq

def find_repeats(seq, k_min=3, k_max=6, min_repeats=2):
    """Find motifs of length k that appear ≥ min_repeats times (overlaps allowed)."""
    n = len(seq)
    results = []
    for k in range(k_min, k_max + 1):
        index = defaultdict(list)
        for i in range(0, n - k + 1):
            kmer = seq[i:i+k]
            index[kmer].append(i)
        for motif, starts in index.items():
            if len(starts) >= min_repeats:
                results.append({
                    "k": k,
                    "motif": motif,
                    "count": len(starts),
                    "positions": [p + 1 for p in starts]
                })
    results.sort(key=lambda d: (d["k"], -d["count"], d["motif"]))
    return results

def print_and_save(header, seq, results, min_repeats, out_file):
    """Print and save summary of results to a text file."""
    with open(out_file, "w") as f:
        f.write(f"{header}\n")
        f.write(f"Sequence length: {len(seq)} bases\n")
        f.write(f"Motifs of length 3–6 appearing ≥ {min_repeats} times\n\n")

        print(f"Header: {header}")
        print(f"Length: {len(seq)} bases\n")
        print(f"Motifs of length 3–6 appearing ≥ {min_repeats} times:\n")

        if not results:
            print("No repeats found.")
            f.write("No repeats found.\n")
            return

        current_k = None
        for r in results:
            if r["k"] != current_k:
                current_k = r["k"]
                print(f"\n{current_k}-base motifs:")
                f.write(f"\n{current_k}-base motifs:\n")
            positions = ", ".join(map(str, r["positions"][:10]))
            if len(r["positions"]) > 10:
                positions += " ..."
            line = f"  {r['motif']} | count={r['count']} | positions: {positions}"
            print(line)
            f.write(line + "\n")

    print(f"\nResults saved to: {out_file}")

def main():
    header, seq = read_fasta(FASTA_FILE)
    if not (1000 <= len(seq) <= 3000):
        print(f"Warning: sequence length {len(seq)} is outside 1000–3000 nt.")
    results = find_repeats(seq, K_MIN, K_MAX, MIN_REPEATS)
    print_and_save(header, seq, results, MIN_REPEATS, OUT_FILE)

if __name__ == "__main__":
    main()

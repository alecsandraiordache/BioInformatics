#download 10 influenza genomes. for each genome plot on chart the most freq repetitions



import os
import re
from collections import defaultdict
import matplotlib.pyplot as plt

FOLDER = "virus_genomes"  
K_MIN, K_MAX = 3, 6       
MIN_REPEATS = 2           


def read_fasta(path):
    with open(path, "r") as f:
        lines = [l.strip() for l in f if l.strip()]
    header = lines[0] if lines and lines[0].startswith(">") else os.path.basename(path)
    seq = "".join(l for l in lines if not l.startswith(">")).upper()
    seq = seq.replace("U", "T")  
    seq = re.sub(r"[^ACGT]", "", seq)  
    return header, seq


def find_repeats(seq, k_min=3, k_max=6, min_repeats=2):
    """Find all motifs (length k) that occur at least N times (overlaps allowed)."""
    n = len(seq)
    results = []
    for k in range(k_min, k_max + 1):
        index = defaultdict(list)
        for i in range(0, n - k + 1):
            kmer = seq[i:i + k]
            index[kmer].append(i)
        for motif, starts in index.items():
            if len(starts) >= min_repeats:
                results.append({
                    "k": k,
                    "motif": motif,
                    "count": len(starts),
                })
    results.sort(key=lambda d: (-d["count"], d["k"], d["motif"]))
    return results


def pick_top_motif(results):
    """Return the single most frequent motif."""
    return results[0] if results else None


def main():
    genomes = []
    for file in sorted(os.listdir(FOLDER)):
        if file.endswith(".fasta"):
            path = os.path.join(FOLDER, file)
            header, seq = read_fasta(path)
            results = find_repeats(seq, K_MIN, K_MAX, MIN_REPEATS)
            top = pick_top_motif(results)
            genomes.append((file, header, len(seq), top))

    print("\n=== Most frequent 3–6 base repeats per genome ===\n")
    for name, header, length, top in genomes:
        if top:
            print(f"{name} ({length} bp): {top['motif']} (k={top['k']}) ×{top['count']}")
        else:
            print(f"{name} ({length} bp): no repeats found")


    labels, values = [], []
    for name, _, _, top in genomes:
        if top:
            labels.append(f"{name}\n{top['motif']}(k={top['k']})")
            values.append(top["count"])
        else:
            labels.append(f"{name}\n(none)")
            values.append(0)

    plt.figure(figsize=(12, 6))
    plt.bar(range(len(values)), values)
    plt.xticks(range(len(values)), labels, rotation=45, ha="right")
    plt.ylabel("Most frequent motif count")
    plt.title("Top 3–6 bp repeats per viral genome")
    plt.tight_layout()
    plt.show()



if __name__ == "__main__":
    main()

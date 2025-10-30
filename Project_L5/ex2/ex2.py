#Search and upload from NCBI on other websites 10 virus viral genomes.
#Use each of these genomes in order to take samples and make a set of samples for each virus
#A. Measure the time (in ms) of assembly for each set
#B. Measure their overall C+G percentage
#C. Plot a chart in which the coord of the points are as follows: 
    #on the Y-axis, the point will take the time in ms
    #on the X-axis, the point will take the overall C+G percentage
#D. Make a text file in which you explain the differences between the positions of the points.


import os
import time
import random
from typing import List, Tuple, Dict
import matplotlib.pyplot as plt


FOLDER = "virus_genomes"   
N_READS = 2000
READ_LEN_RANGE: Tuple[int, int] = (100, 150)
MIN_OVERLAP = 10
SEED = 123
random.seed(SEED)


def read_fasta(path: str) -> str:
    seq_parts = []
    with open(path, "r") as f:
        for line in f:
            if not line.startswith(">"):
                seq_parts.append(line.strip().upper())
    return "".join(seq_parts)

def sample_reads(genome: str, n_reads: int, len_range: Tuple[int, int]) -> list:
    L = len(genome)
    reads = []
    for _ in range(n_reads):
        rlen = random.randint(len_range[0], len_range[1])
        start = random.randint(0, L - rlen)
        reads.append(genome[start:start+rlen])
    return reads

def overlap_len(a: str, b: str, min_len: int) -> int:
    seed = b[:min_len]
    start = 0
    while True:
        pos = a.find(seed, start)
        if pos == -1:
            return 0
        if b.startswith(a[pos:]):
            return len(a) - pos
        start = pos + 1

def build_prefix_index(reads: list, k: int) -> dict:
    idx = {}
    for i, s in enumerate(reads):
        if len(s) >= k:
            key = s[:k]
            if key not in idx:
                idx[key] = []
            idx[key].append(i)
    return idx

def assemble_greedy(reads: list, min_overlap: int) -> str:
    n = len(reads)
    used = [False] * n
    k = min_overlap
    idx = build_prefix_index(reads, k)
    contigs = []

    for i in range(n):
        if used[i]:
            continue
        contig = reads[i]
        used[i] = True
        while len(contig) >= k:
            suffix = contig[-k:]
            cand = idx.get(suffix, [])
            best_j, best_ol = -1, 0
            for j in cand:
                if used[j] or j == i:
                    continue
                ol = overlap_len(contig, reads[j], k)
                if ol > best_ol:
                    best_ol, best_j = ol, j
            if best_ol >= k and best_j != -1:
                contig += reads[best_j][best_ol:]
                used[best_j] = True
            else:
                break
        contigs.append(contig)

    return max(contigs, key=len) if contigs else ""

def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    return 100 * (seq.count("G") + seq.count("C")) / len(seq)


def main():
    viruses = [f for f in os.listdir(FOLDER) if f.endswith(".fasta")]
    results = []

    for f in viruses:
        name = os.path.splitext(f)[0]
        path = os.path.join(FOLDER, f)
        print(f"[{name}] Reading genome...")
        genome = read_fasta(path)
        reads = sample_reads(genome, N_READS, READ_LEN_RANGE)

        start = time.perf_counter()
        contig = assemble_greedy(reads, MIN_OVERLAP)
        elapsed = (time.perf_counter() - start) * 1000

        gc = gc_percent(genome)
        results.append((name, gc, elapsed))
        print(f"  GC%={gc:.2f}  Time={elapsed:.1f} ms  Length={len(genome)} bp")

    xs = [r[1] for r in results]
    ys = [r[2] for r in results]
    labels = [r[0] for r in results]

    plt.figure(figsize=(8, 6))
    plt.scatter(xs, ys)
    for x, y, label in zip(xs, ys, labels):
        plt.annotate(label, (x, y), textcoords="offset points", xytext=(5, 5), fontsize=9)
    plt.xlabel("Overall C+G (%)")
    plt.ylabel("Assembly time (ms)")
    plt.title("Assembly time vs GC% (10 viral genomes)")
    plt.tight_layout()
    plt.show()

    with open("notes.txt", "w") as f:
        f.write("Assembly Time vs GC% — Notes\n")
        f.write("--------------------------------\n\n")
        for n, gc, t in results:
            f.write(f"{n:15s}  GC%={gc:.2f}  Time={t:.1f} ms\n")
        f.write("\nDifferences explanation:\n")
        f.write("Viruses with extreme GC% or longer genomes may take longer assembly times.\n")
        f.write("Low coverage (relative to genome length) or repetitive regions increase overlap checks.\n")
        f.write("GC-rich genomes have higher sequence diversity → faster overlap matching.\n")
        f.write("AT-rich or repetitive genomes create ambiguous overlaps → slower assembly.\n")

    print("\n[✔] Done! Results saved in notes.txt")

if __name__ == "__main__":
    main()

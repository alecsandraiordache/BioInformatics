#1. Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology),
#between 1000 and 3000 nucleotides (letters).

#2. Take 2000 random samples from this sequence, of about 100-150 bases.

#3. Store these samples in an array.

#4. Rebuild the original DNA sequence using these random samples.

#What would be the main problem with the algorithm approach ? Describe; note what kind of struct 
#inside the original may create different computation oscos

#note : the samples must be aligned starting with the min of 10 possitions in order to
# avoid random matching

import os
import random
from typing import List, Tuple, Dict

FASTA_FILE = "sequence.fasta"
N_READS = 2000
READ_LEN_RANGE: Tuple[int, int] = (100, 150)
MIN_OVERLAP = 10
SEED = 42
random.seed(SEED)

def read_fasta(path: str) -> str:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing '{path}'. Put the FASTA in the same folder as this script.")
    seq_parts: List[str] = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith(">"):
                continue
            seq_parts.append(line.strip().upper())
    seq = "".join(seq_parts)
    if not seq:
        raise ValueError("FASTA has no sequence lines (only headers?)")
    if len(seq) < READ_LEN_RANGE[1]:
        raise ValueError(f"Sequence too short ({len(seq)} bp). Needs at least {READ_LEN_RANGE[1]} bp.")
    return seq

def sample_reads(genome: str, n_reads: int, len_range: Tuple[int, int]) -> List[str]:
    L = len(genome)
    reads: List[str] = []
    for _ in range(n_reads):
        rlen = random.randint(len_range[0], len_range[1])
        start = random.randint(0, L - rlen)
        reads.append(genome[start:start + rlen])
    return reads

def overlap_len(a: str, b: str, min_len: int) -> int:
    """Longest suffix of a that matches prefix of b (>= min_len)."""
    seed = b[:min_len]
    start = 0
    while True:
        pos = a.find(seed, start)
        if pos == -1:
            return 0
        if b.startswith(a[pos:]):
            return len(a) - pos
        start = pos + 1

def build_prefix_index(reads: List[str], k: int) -> Dict[str, List[int]]:
    idx: Dict[str, List[int]] = {}
    for i, s in enumerate(reads):
        if len(s) >= k:
            idx.setdefault(s[:k], []).append(i)
    return idx

def assemble_greedy_kmer(reads: List[str], min_overlap: int) -> str:
    """Greedy assembly using a k-mer (k=min_overlap) prefix index; returns the longest contig."""
    n = len(reads)
    used = [False] * n
    k = min_overlap
    prefix_idx = build_prefix_index(reads, k)
    contigs: List[str] = []

    for i in range(n):
        if used[i]:
            continue
        contig = reads[i]
        used[i] = True
        while len(contig) >= k:
            suffix = contig[-k:]
            candidates = prefix_idx.get(suffix, [])
            best_j = -1
            best_ol = 0
            for j in candidates:
                if used[j] or j == i:
                    continue
                ol = overlap_len(contig, reads[j], k)
                if ol > best_ol:
                    best_ol = ol
                    best_j = j
            if best_ol >= k and best_j != -1:
                contig += reads[best_j][best_ol:]
                used[best_j] = True
            else:
                break
        contigs.append(contig)

    return max(contigs, key=len) if contigs else ""

def main():
    # 1)
    original_sequence = read_fasta(FASTA_FILE)

    # 2)
    samples = sample_reads(original_sequence, N_READS, READ_LEN_RANGE)


    # 4) 
    reconstructed_sequence = assemble_greedy_kmer(samples, MIN_OVERLAP)

    print(" DONE ")
    print(f"Original length     : {len(original_sequence)} bp")
    print(f"Reconstructed length: {len(reconstructed_sequence)} bp")

if __name__ == "__main__":
    main()

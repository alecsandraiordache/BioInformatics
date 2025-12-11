#download from NCBI the influenza genomes and covid 19 genomes and align these genomes by using 
#the local alignment method 
#Note that you have to add an in between layer solution in order to be able to align 
#the 2 genomes files. 
# The alignment made step by step on big region with the connection between the results
# Note that the seq align. alg. doesn t allow for the aligh. of seqs. because the bigger seq, 
#the bigger scoring matrix
#The main result should be thevizualization of the similarities between the 2 genomes
#NOTE: Pls do not use the shortcuts given by AI . Use native code



import os
import numpy as np
import matplotlib.pyplot as plt



WINDOW_SIZE = 300  
STEP = 300  

MATCH = 2
MISMATCH = -1
GAP = -2



def read_first_fasta(folder):
    files = sorted([
        f for f in os.listdir(folder)
        if f.lower().endswith(".fasta") or f.lower().endswith(".fa")
    ])
    path = os.path.join(folder, files[0])
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"): 
                continue
            seq.append(line.upper())
    sequence = "".join(seq)
    print(f"Loaded {len(sequence)} bases from {path}")
    return sequence



def smith_waterman_score_only(s1, s2):
    n, m = len(s1), len(s2)
    H = [[0]*(m+1) for _ in range(n+1)]
    best = 0

    for i in range(1, n+1):
        ch1 = s1[i-1]
        row = H[i]
        prev = H[i-1]
        for j in range(1, m+1):
            ch2 = s2[j-1]

            diag = prev[j-1] + (MATCH if ch1 == ch2 else MISMATCH)
            up   = prev[j] + GAP
            left = row[j-1] + GAP

            v = max(0, diag, up, left)
            row[j] = v
            if v > best:
                best = v
    return best



def windows(seq, size, step):
    pts = []
    start = 0
    n = len(seq)
    while start < n:
        end = min(n, start + size)
        pts.append((start, end))
        if end == n:
            break
        start += step
    return pts



def build_matrix(covid, flu):
    covid_win = windows(covid, WINDOW_SIZE, STEP)
    flu_win   = windows(flu, WINDOW_SIZE, STEP)

    C = len(covid_win)
    F = len(flu_win)

    M = np.zeros((C, F), dtype=int)

    for i, (cs, ce) in enumerate(covid_win):
        frag_c = covid[cs:ce]
        for j, (fs, fe) in enumerate(flu_win):
            frag_f = flu[fs:fe]
            M[i, j] = smith_waterman_score_only(frag_c, frag_f)

    return M



def plot_map(matrix):
    plt.figure(figsize=(10,6))
    plt.imshow(matrix, cmap="inferno", aspect="auto", origin="lower")
    plt.colorbar(label="Smith–Waterman Score")
    plt.xlabel("Influenza Windows")
    plt.ylabel("COVID Windows")
    plt.title("COVID vs Influenza — Local Alignment Similarity Map")
    plt.tight_layout()
    plt.show()



def main():
    covid = read_first_fasta("covid")
    flu   = read_first_fasta("influenza")

    M = build_matrix(covid, flu)
    plot_map(M)


if __name__ == "__main__":
    main()

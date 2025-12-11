#Formulate 3 scorring eq. that are able to show the level of similarity between the 2 seq aligned
#in the previous asignment 
#Implement each of these 3 scoring eqs. in your current implementation 



import os
import numpy as np
import matplotlib.pyplot as plt

WINDOW_SIZE = 300
STEP = 300

MATCH = 2
MISMATCH = -1
GAP = -2

def read_first_fasta(folder):
    files = sorted([f for f in os.listdir(folder)
                    if f.lower().endswith(".fasta") or f.lower().endswith(".fa")])
    if not files:
        raise FileNotFoundError(f"No FASTA files in {folder}")

    path = os.path.join(folder, files[0])
    seq = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"): continue
            seq.append(line.upper())

    sequence = "".join(seq)
    print(f"Loaded {len(sequence)} bases from {path}")
    return sequence



def smith_waterman_metrics(s1, s2):
    n, m = len(s1), len(s2)
    H  = [[0]*(m+1) for _ in range(n+1)]
    TB = [[0]*(m+1) for _ in range(n+1)]

    max_score = 0
    max_i = 0
    max_j = 0

    for i in range(1, n+1):
        row  = H[i]
        prow = H[i-1]
        tb_row = TB[i]
        ch1 = s1[i-1]

        for j in range(1, m+1):
            ch2 = s2[j-1]

            diag = prow[j-1] + (MATCH if ch1 == ch2 else MISMATCH)
            up   = prow[j] + GAP
            left = row[j-1] + GAP

            v = max(0, diag, up, left)

            if v == 0:
                tb = 0
            elif v == diag:
                tb = 1
            elif v == up:
                tb = 2
            else:
                tb = 3

            row[j] = v
            tb_row[j] = tb

            if v > max_score:
                max_score = v
                max_i = i
                max_j = j

    i, j = max_i, max_j
    M = X = G = L = 0

    while i > 0 and j > 0 and TB[i][j] != 0:
        step = TB[i][j]

        if step == 1: 
            L += 1
            if s1[i-1] == s2[j-1]:
                M += 1
            else:
                X += 1
            i -= 1
            j -= 1

        elif step == 2: 
            G += 1
            L += 1
            i -= 1

        elif step == 3: 
            G += 1
            L += 1
            j -= 1

    if L == 0:
        return max_score, 0.0, 0.0, 0.0

    PID = M / L

    CAI = M / min(n, m)

    NSW = max_score / (2 * L)

    return max_score, PID, CAI, NSW



def windows(seq, size, step):
    result = []
    start = 0
    n = len(seq)
    while start < n:
        end = min(n, start + size)
        result.append((start, end))
        if end == n:
            break
        start += step
    return result



def build_matrices(covid, flu):
    covid_win = windows(covid, WINDOW_SIZE, STEP)
    flu_win   = windows(flu, WINDOW_SIZE, STEP)

    C = len(covid_win)
    F = len(flu_win)

    print(f"COVID windows     : {C}")
    print(f"Influenza windows : {F}")

    M_pid = np.zeros((C, F), dtype=float)
    M_cai = np.zeros((C, F), dtype=float)
    M_nsw = np.zeros((C, F), dtype=float)

    for i, (cs, ce) in enumerate(covid_win):
        frag_c = covid[cs:ce]
        print(f"Processing window {i+1}/{C}...")

        for j, (fs, fe) in enumerate(flu_win):
            frag_f = flu[fs:fe]
            _, pid, cai, nsw = smith_waterman_metrics(frag_c, frag_f)

            M_pid[i, j] = pid
            M_cai[i, j] = cai
            M_nsw[i, j] = nsw

    return M_pid, M_cai, M_nsw



def plot_maps(M_pid, M_cai, M_nsw):

    fig, axes = plt.subplots(1, 3, figsize=(17, 5), sharey=True)

    im0 = axes[0].imshow(M_pid, cmap="viridis", aspect="auto", origin="lower")
    axes[0].set_title("PID (Percent Identity)")
    axes[0].set_xlabel("Influenza windows")
    axes[0].set_ylabel("COVID windows")
    fig.colorbar(im0, ax=axes[0], fraction=0.046)

    im1 = axes[1].imshow(M_cai, cmap="plasma", aspect="auto", origin="lower")
    axes[1].set_title("CAI (Coverage Adjusted Identity)")
    axes[1].set_xlabel("Influenza windows")
    fig.colorbar(im1, ax=axes[1], fraction=0.046)

    im2 = axes[2].imshow(M_nsw, cmap="magma", aspect="auto", origin="lower")
    axes[2].set_title("NSW (Normalized SW score)")
    axes[2].set_xlabel("Influenza windows")
    fig.colorbar(im2, ax=axes[2], fraction=0.046)

    fig.suptitle("COVID vs Influenza — Similarity Maps (3 Different Scoring Equations)")
    plt.tight_layout()
    plt.show()


def main():
    covid = read_first_fasta("covid")
    flu   = read_first_fasta("influenza")

    M_pid, M_cai, M_nsw = build_matrices(covid, flu)
    plot_maps(M_pid, M_cai, M_nsw)


if __name__ == "__main__":
    main()

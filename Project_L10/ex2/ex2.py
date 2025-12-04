# Download 10 influenza genomes and 10 covid-19 genomes. Use their FASTA FILES 
#to plor rheir objective digital staints. On a 2nd type of chart, plot the center of weigth for 
#each ODS. Lable each of the centers of weight with the name and variant of the virus

import os
import glob
import statistics
import matplotlib.pyplot as plt

def read_fasta_single(path: str) -> str:
    seq_lines = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_lines.append(line)
    return "".join(seq_lines).upper()


def sliding_windows(seq: str, win_size: int):
    L = len(seq)
    return [seq[i:i+win_size] for i in range(0, L - win_size + 1)]

def cg_percent(seq: str) -> float:
    if not seq:
        return 0.0
    cg = sum(1 for b in seq if b in ("C", "G"))
    return (cg / len(seq)) * 100.0

def kappa_ic(seq: str) -> float:
    L = len(seq)
    if L < 2:
        return 0.0

    N = L - 1
    total = 0.0

    for u in range(1, N + 1):
        matches = 0
        B_len = L - u
        for i in range(B_len):
            if seq[i] == seq[u + i]:
                matches += 1
        total += (matches / B_len) * 100.0

    return total / N

def compute_pattern(seq: str, win_size: int):
    windows = sliding_windows(seq, win_size)
    xs = [cg_percent(w) for w in windows]
    ys = [kappa_ic(w) for w in windows]
    return xs, ys

def center_of_weight(xs, ys):
    if not xs or not ys:
        return 0.0, 0.0
    return statistics.mean(xs), statistics.mean(ys)


def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))

    covid_dir = os.path.join(base_dir, "covid")
    influenza_dir = os.path.join(base_dir, "influenza")

    covid_files = sorted(glob.glob(os.path.join(covid_dir, "*.fasta")))
    influenza_files = sorted(glob.glob(os.path.join(influenza_dir, "*.fasta")))

    print("Found COVID files:")
    for f in covid_files:
        print("  ", os.path.basename(f))

    print("\nFound influenza files:")
    for f in influenza_files:
        print("  ", os.path.basename(f))

    win_size = 30

    all_centers = []
    all_labels = []
    all_types = []  
    
    plt.figure()
    plt.title("Objective Digital Stains - Influenza vs COVID-19")
    plt.xlabel("C+G%")
    plt.ylabel("Kappa IC")
    plt.grid(True)

    for idx, path in enumerate(covid_files, start=1):
        seq = read_fasta_single(path)
        xs, ys = compute_pattern(seq, win_size)
        cx, cy = center_of_weight(xs, ys)

        all_centers.append((cx, cy))
        label = f"COVID-{idx}"
        all_labels.append(label)
        all_types.append("covid")

        plt.scatter(xs, ys, s=5, alpha=0.5, marker="s", label=None)

        print(f"COVID-{idx}: length={len(seq)}, center=({cx:.2f}, {cy:.2f})")

    for idx, path in enumerate(influenza_files, start=1):
        seq = read_fasta_single(path)
        xs, ys = compute_pattern(seq, win_size)
        cx, cy = center_of_weight(xs, ys)

        all_centers.append((cx, cy))
        label = f"Influenza-{idx}"
        all_labels.append(label)
        all_types.append("influenza")

        plt.scatter(xs, ys, s=5, alpha=0.5, marker="o", label=None)

        print(f"Influenza-{idx}: length={len(seq)}, center=({cx:.2f}, {cy:.2f})")

    plt.figure()
    plt.title("Centers of Weight of ODS")
    plt.xlabel("C+G% (center)")
    plt.ylabel("Kappa IC (center)")
    plt.grid(True)

    for (cx, cy), label, vtype in zip(all_centers, all_labels, all_types):
        if vtype == "covid":
            plt.scatter(cx, cy, marker="s")
        else:
            plt.scatter(cx, cy, marker="o")
        plt.text(cx, cy, "  " + label, fontsize=8)

    plt.show()

if __name__ == "__main__":
    main()

from pathlib import Path
import re
import math
import matplotlib.pyplot as plt


FASTA_FILE = "pUC19.fasta"
IS_CIRCULAR = True   

enzymes = {
    "EcoRI": {
        "site": "GAATTC",
        "cut_offset": 1      # G|AATTC
    },
    "BamHI": {
        "site": "GGATCC",
        "cut_offset": 1      # G|GATCC
    },
    "HindIII": {
        "site": "AAGCTT",
        "cut_offset": 1      # A|AGCTT
    },
    "TaqI": {
        "site": "TCGA",
        "cut_offset": 1      # T|CGA
    },
    "HaeIII": {
        "site": "GGCC",
        "cut_offset": 2      # GG|CC
    },
}


def read_fasta_sequence(path: str) -> str:
    text = Path(path).read_text().strip().splitlines()
    if not text or not text[0].startswith(">"):
        raise ValueError("Invalid FASTA file.")

    header = text[0]
    seq = "".join(line.strip().upper() for line in text[1:]
                  if not line.startswith(">"))

    for ch in seq:
        if ch not in "ACGTN":
            raise ValueError(f"Invalid nucleotide found: {ch}")

    print(f"FASTA file loaded: {path}")
    print(f"Header: {header}")
    print(f"Sequence length: {len(seq)} bp\n")
    return seq


def find_cut_positions(seq: str, motif: str, cut_offset: int):
    pattern = f"(?=({motif}))"  # regex lookahead
    positions = []
    for m in re.finditer(pattern, seq):
        start = m.start()
        cut_pos = start + cut_offset
        positions.append(cut_pos)
    return positions


def compute_fragments(seq_length: int, cuts: list, circular: bool):
    if not cuts:
        return [seq_length]

    cuts = sorted(cuts)

    if circular and len(cuts) == 1:
        return [seq_length]

    if not circular:
        fragments = []
        last = 0
        for c in cuts:
            fragments.append(c - last)
            last = c
        fragments.append(seq_length - last)
        return fragments

    fragments = []
    n = len(cuts)
    for i in range(n):
        start = cuts[i]
        end = cuts[(i + 1) % n]
        if end >= start:
            fragments.append(end - start)
        else:
            fragments.append(seq_length - start + end)
    return fragments


def print_digest_report(enzyme_name: str, seq_len: int, cuts: list, fragments: list):
    print(f"{enzyme_name}")
    print(f"Sequence length: {seq_len} bp")
    print(f"Number of cuts: {len(cuts)}")

    if cuts:
        print("Cut positions (1-based):", [c + 1 for c in sorted(cuts)])
    else:
        print("No recognition sites found.")

    print("Fragment sizes (bp):", sorted(fragments, reverse=True))
    print()



def plot_gel(fragments_by_enzyme: dict):
    sizes = [s for frags in fragments_by_enzyme.values() for s in frags if s > 0]
    if not sizes:
        print("No fragments to plot.")
        return

    max_size = max(sizes)
    min_size = min(sizes)
    max_log = math.log10(max_size)
    min_log = math.log10(min_size)

    fig, ax = plt.subplots(figsize=(6, 6))

    enzyme_names = list(fragments_by_enzyme.keys())
    lane_positions = range(len(enzyme_names))  

    for lane, enz in zip(lane_positions, enzyme_names):
        frags = fragments_by_enzyme[enz]
        if not frags:
            continue

        for s in sorted(frags):
            if s <= 0:
                continue

            log_s = math.log10(s)
            frac = (max_log - log_s) / (max_log - min_log + 1e-9)
            y = frac

            x_left = lane - 0.35
            x_right = lane + 0.35
            ax.hlines(y, x_left, x_right, linewidth=6)

    ax.set_xticks(list(lane_positions))
    ax.set_xticklabels(enzyme_names)

    ax.set_ylabel("Migration distance (smaller fragments migrate further)")
    ax.set_ylim(1.05, -0.05)

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    plt.title("Restriction Digest – Electrophoresis Gel Simulation")
    plt.tight_layout()
    plt.show()


def main():
    seq = read_fasta_sequence(FASTA_FILE)
    seq_len = len(seq)

    fragments_by_enzyme = {}

    for enz_name, info in enzymes.items():
        site = info["site"]
        offset = info["cut_offset"]

        cuts = find_cut_positions(seq, site, offset)
        fragments = compute_fragments(seq_len, cuts, IS_CIRCULAR)

        fragments_by_enzyme[enz_name] = fragments
        print_digest_report(enz_name, seq_len, cuts, fragments)

    plot_gel(fragments_by_enzyme)


if __name__ == "__main__":
    main()

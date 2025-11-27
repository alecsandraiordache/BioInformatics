#Download 10 influenza viruses variants from NCBI and analyze their genome by using the 
#application from assigment 1.
#a.  Make an electrophoresys for each genome
#b. eliminate all lines that are in common between the gel simullation such that the 
#differences will be show.
#Merge all the electrophoresys gel simulation (that show ONLY the differences) in one general el. gel.


from pathlib import Path
import re
import math
import matplotlib.pyplot as plt
import sys


LOGFILE = "output.txt"

class Tee(object):
    def __init__(self, filename):
        self.file = open(filename, "w")
        self.stdout = sys.stdout

    def write(self, message):
        self.stdout.write(message)
        self.file.write(message)

    def flush(self):
        self.stdout.flush()
        self.file.flush()


sys.stdout = Tee(LOGFILE)


GENOME_FOLDER = Path("virus_genomes")
IS_CIRCULAR = False
COMMON_TOLERANCE = 10

enzymes = {
    "EcoRI":  {"site": "GAATTC", "cut_offset": 1},
    "BamHI":  {"site": "GGATCC", "cut_offset": 1},
    "HindIII":{"site": "AAGCTT", "cut_offset": 1},
    "TaqI":   {"site": "TCGA",   "cut_offset": 1},
    "HaeIII": {"site": "GGCC",   "cut_offset": 2},
}


def read_fasta_sequence(path: Path) -> str:
    allowed = set("ACGTURYSWKMBDHVN")  # full IUPAC

    lines = path.read_text().strip().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"Invalid FASTA file: {path}")

    seq = "".join(
        line.strip().upper()
        for line in lines[1:]
        if not line.startswith(">")
    )

    for ch in seq:
        if ch not in allowed:
            raise ValueError(f"Invalid nucleotide '{ch}' in file: {path.name}")

    return seq.replace("U", "T")


def find_cut_positions(seq, motif, cut_offset):
    pattern = f"(?=({motif}))"
    return [m.start() + cut_offset for m in re.finditer(pattern, seq)]

def compute_fragments(length, cuts, circular):
    if not cuts:
        return [length]

    cuts = sorted(cuts)

    if circular and len(cuts) == 1:
        return [length]

    if not circular:
        frags, last = [], 0
        for c in cuts:
            frags.append(c - last)
            last = c
        frags.append(length - last)
        return frags

    frags = []
    for i in range(len(cuts)):
        start = cuts[i]
        end = cuts[(i+1) % len(cuts)]
        if end >= start:
            frags.append(end - start)
        else:
            frags.append(length - start + end)
    return frags

def digest_genome(seq):
    length = len(seq)
    results = {}

    for enz, info in enzymes.items():
        cuts = find_cut_positions(seq, info["site"], info["cut_offset"])
        frags = compute_fragments(length, cuts, IS_CIRCULAR)
        results[enz] = frags

        print(f"Genome length: {length} bp | Enzyme: {enz}")
        print(f"  Cuts: {len(cuts)}")
        print(f"  Fragments: {sorted(frags, reverse=True)}")

    print()
    return results


def is_close(a, b, tol): return abs(a-b) <= tol

def find_common(fragment_lists, tol):
    base = sorted(fragment_lists[0])
    common = []
    for x in base:
        if all(any(is_close(x, y, tol) for y in lst) for lst in fragment_lists[1:]):
            common.append(x)
    return common

def remove_common(fragment_lists, common, tol):
    out = []
    for lst in fragment_lists:
        filtered = [x for x in lst if not any(is_close(x, c, tol) for c in common)]
        out.append(filtered)
    return out


def plot_gel(enz, filtered_lists, labels):
    sizes = [s for lst in filtered_lists for s in lst]
    if not sizes:
        print(f"No differing bands for enzyme {enz}.")
        return

    max_s, min_s = max(sizes), min(sizes)
    max_log, min_log = math.log10(max_s), math.log10(min_s)

    fig, ax = plt.subplots(figsize=(6, 8))

    for lane, (label, frags) in enumerate(zip(labels, filtered_lists)):
        for s in sorted(frags):
            log_s = math.log10(s)
            frac = (max_log - log_s) / (max_log - min_log + 1e-9)
            y = 0.05 + frac * 0.9
            ax.hlines(y, lane - 0.35, lane + 0.35, linewidth=8)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticks([])
    ax.set_ylim(1, 0)

    for side in ["top", "right", "left"]:
        ax.spines[side].set_visible(False)

    ax.set_title(f"Gel (differences only) – {enz}")
    plt.tight_layout()
    plt.show()


def main():
    fasta_files = sorted(GENOME_FOLDER.glob("*.fasta"))
    if not fasta_files:
        raise SystemExit("No FASTA files in virus_genomes/")

    print("Virus genomes detected")
    for f in fasta_files:
        print("  ", f.name)
    print()

    all_results = []
    labels = []

    for fasta in fasta_files:
        print(f"Processing {fasta.name} ")
        seq = read_fasta_sequence(fasta)
        labels.append(fasta.stem)
        res = digest_genome(seq)
        all_results.append(res)

    for enz in enzymes:
        print(f"Analyzing enzyme {enz}")
        frag_lists = [r[enz] for r in all_results]
        common = find_common(frag_lists, COMMON_TOLERANCE)
        print("Common bands:", sorted(common))
        filtered = remove_common(frag_lists, common, COMMON_TOLERANCE)
        plot_gel(enz, filtered, labels)

if __name__ == "__main__":
    main()

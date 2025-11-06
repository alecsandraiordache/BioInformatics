#Download 10 influenza virus genomes and apply the electrophoresis gel simulation an each of them
# Make a comparison between the 10 electroph. gel simulations and show which of the influenza 
#genomes show the most DNA segments. You can plot them in the same graph but also separately 
#because may overlap. 
#As the main restriction enzimes, please use ECOR1

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def read_fasta(path):
    with open(path, "r") as f:
        seq = "".join([ln.strip() for ln in f if not ln.startswith(">")]).upper()
    return "".join([c for c in seq if c in ("A", "C", "G", "T", "N")])


def digest_ecori(seq):
    motif = "GAATTC"
    cut_offset = 1
    cut_sites = []
    i = 0
    while True:
        j = seq.find(motif, i)
        if j == -1:
            break
        cut_sites.append(j + cut_offset)
        i = j + 1
    if not cut_sites:
        return [len(seq)]
    positions = [0] + cut_sites + [len(seq)]
    fragments = [positions[k+1] - positions[k] for k in range(len(positions)-1)]
    return [f for f in fragments if f > 0]


def migration_distance(bp, bp_min, bp_max):
    log_min = np.log10(bp_min)
    log_max = np.log10(bp_max)
    return (np.log10(bp_max) - np.log10(bp)) / (log_max - log_min + 1e-9)


folder = "virus_genomes"
paths = sorted(glob.glob(os.path.join(folder, "*.fasta")))
if not paths:
    raise FileNotFoundError("No .fasta files found in ./virus_genomes")

genome_names = [os.path.splitext(os.path.basename(p))[0] for p in paths]
sequences = [read_fasta(p) for p in paths]
fragments = [digest_ecori(seq) for seq in sequences]


band_counts = [len(f) for f in fragments]
max_count = max(band_counts)
most_fragments = [genome_names[i] for i, c in enumerate(band_counts) if c == max_count]

print("=== EcoRI fragments per genome ===")
for name, count in zip(genome_names, band_counts):
    print(f"{name}: {count} fragments")

print("\nGenome(s) with the most fragments:")
for name in most_fragments:
    print(f"- {name}")


all_lengths = [bp for lst in fragments for bp in lst]
bp_min = max(50, min(all_lengths))
bp_max = max(all_lengths)

fig, ax = plt.subplots(figsize=(1.5 * len(paths), 8))

ladder_sizes = [100, 200, 300, 500, 700, 1000, 1500, 2000, 3000]
for b in ladder_sizes:
    if bp_min <= b <= bp_max:
        y = 0.08 + 0.86 * migration_distance(b, bp_min, bp_max)
        ax.plot([0.5, 0.8], [y, y], linewidth=3)
        ax.text(0.45, y, f"{b}", ha="right", va="center", fontsize=8)


for i, (name, bands) in enumerate(zip(genome_names, fragments), start=1):
    x_center = i + 1 
    ax.add_patch(Rectangle((x_center - 0.3, 0.05), 0.6, 0.90, alpha=0.2))
    ax.add_patch(Rectangle((x_center - 0.28, 0.94), 0.56, 0.02, fill=False, linewidth=1.0))
    for bp in bands:
        y = 0.08 + 0.86 * migration_distance(bp, bp_min, bp_max)
        ax.plot([x_center - 0.22, x_center + 0.22], [y, y], linewidth=4)
    ax.text(x_center, 0.01, name, ha="center", va="bottom", rotation=90, fontsize=8)

ax.set_xlim(0, len(paths) + 2)
ax.set_ylim(0, 1)
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("EcoRI Digest – Influenza Genomes (All Lanes Combined)")
plt.tight_layout()
plt.show()

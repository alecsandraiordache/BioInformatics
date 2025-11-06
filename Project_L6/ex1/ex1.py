#Gel electrophoresis is an analysis method implemented in all disciplines of life sciences. 
#The results of gel electrophoresis indicate the relative sizes of fragments, 
#which is useful for restriction mapping and analyzing PCR fragments. 

#1. Take an arbitrary DNA sequence from the NCBI (National Center for Biotechnology), between 1000 and 3000 nucleotides (letters).

#2. Take 10 random samples from this sequence, between 100-3000 bases.

#3. Store these samples in an array.

#4. Simulate the migration of these DNA segments on the electrophoresis gel, based on their 
#molecular weights - however, their length should be sufficient for this exercise 
#(show a visual representation).

#Note: Short DNA fragments meet small friction forces and travel faster through the 
#electrophoresis gel. Long DNA fragments exhibit a high friction force and travel slowly
# through the electrophoresis gel.


import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path


fasta_path = Path("covid.fasta") 
with open(fasta_path) as f:
    lines = f.readlines()

dna_seq = "".join([line.strip() for line in lines if not line.startswith(">")]).upper()
print(f"COVID genome length: {len(dna_seq)} bp")


random.seed(42)
num_samples = 10
min_len = 100
max_len = 3000

samples, lengths = [], []

for _ in range(num_samples):
    L = random.randint(min_len, max_len)
    start = random.randint(0, len(dna_seq) - L)
    frag = dna_seq[start:start + L]
    samples.append(frag)
    lengths.append(L)

def norm_migration(bp, min_bp, max_bp):
    log_min = np.log10(min_bp)
    log_max = np.log10(max_bp)
    return (np.log10(max_bp) - np.log10(bp)) / (log_max - log_min + 1e-9)

min_bp = max(50, min(lengths))
max_bp = max(lengths)
mig = [norm_migration(b, min_bp, max_bp) for b in lengths]


fig, ax = plt.subplots(figsize=(4, 8))

ax.add_patch(Rectangle((0.35, 0.05), 0.30, 0.90, facecolor="#f5e6c8", alpha=0.5))

ax.add_patch(Rectangle((0.36, 0.94), 0.28, 0.02, fill=False, linewidth=1.5, edgecolor="black"))


for i, d in enumerate(mig):
    y = 0.08 + 0.86 * d
    ax.plot([0.38, 0.62], [y, y], linewidth=6)
    ax.text(0.64, y, f"{lengths[i]} bp", va="center", fontsize=8)


ladder_sizes = [3000, 2000, 1500, 1000, 700, 500, 300, 200, 100]
ladder_in = [s for s in ladder_sizes if min_bp <= s <= max_bp]
for b in ladder_in:
    d = norm_migration(b, min_bp, max_bp)
    y = 0.08 + 0.86 * d
    ax.plot([0.15, 0.31], [y, y], linewidth=3)
    ax.text(0.14, y, f"{b}", va="center", ha="right", fontsize=8)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("Simulated Gel Electrophoresis – SARS-CoV-2 (covid.fasta)")

plt.show()

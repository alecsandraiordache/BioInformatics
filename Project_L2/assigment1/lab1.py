import itertools
from collections import Counter


S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA".strip().upper()
n = len(S)

nucleotides = ['A', 'C', 'G', 'T']
dinucleotides = [''.join(p) for p in itertools.product(nucleotides, repeat=2)]
trinucleotides = [''.join(p) for p in itertools.product(nucleotides, repeat=3)]

def count_kmers(seq, k):
    """Numără k-merii suprapuși dintr-o secvență."""
    return Counter(seq[i:i + k] for i in range(len(seq) - k + 1))

dinucleotide_counts = count_kmers(S, 2)
trinucleotide_counts = count_kmers(S, 3)

total_dinucleotide = n - 1
total_trinucleotide = n - 2

dinucleotide_percent = {k: dinucleotide_counts.get(k, 0) / total_dinucleotide * 100 for k in dinucleotides}
trinucleotide_percent = {k: trinucleotide_counts.get(k, 0) / total_trinucleotide * 100 for k in trinucleotides}

print(f"Sequence length: {n}")
print(f"Total dinucleotide windows: {total_dinucleotide}")
print(f"Total trinucleotide windows: {total_trinucleotide}\n")

print("Dinucleotide Percentages:")
print("Dinucleotide | Count | Percentage (%)")
for k in dinucleotides:
    print(f"{k:<12} | {dinucleotide_counts.get(k,0):>5} | {dinucleotide_percent[k]:>10.2f}")

print("Trinucleotide Percentages:")
print("Trinucleotide | Count | Percentage (%)")
for k in trinucleotides:
    print(f"{k:<13} | {trinucleotide_counts.get(k,0):>5} | {trinucleotide_percent[k]:>10.2f}")

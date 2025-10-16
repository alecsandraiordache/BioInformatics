import math

# Read DNA sequence from user
seq = input("Enter DNA sequence (A/C/G/T): ").strip().upper().replace(" ", "")

# Validate the sequence
if not seq or any(ch not in "ACGT" for ch in seq):
    print("Error: The sequence must contain only the letters A, C, G, and T.")
else:
    L = len(seq)
    A = seq.count("A")
    T = seq.count("T")
    G = seq.count("G")
    C = seq.count("C")

    # Formula 1: Simple formula
    tm_simple = 4 * (G + C) + 2 * (A + T)

    # Formula 2: Salt-adjusted formula (Na+ = 0.05 M)
    na = 0.05
    gc_percent = 100.0 * (G + C) / L
    tm_salt = 81.5 + 16.6 * math.log10(na) + 0.41 * gc_percent - 600.0 / L

    # Display the results
    print(f"\nTm (simple formula): {tm_simple:.2f} °C")
    print(f"Tm (salt-adjusted formula): {tm_salt:.2f} °C")

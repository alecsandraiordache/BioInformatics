#Download from NCBA the fasta files containing the COVID 19 genome and the influenza genome
#Use the AI to compare the codon freq between the two genomes.
#a. make a chart that shows the top 10 most frew codons for COVID 19
#b. make a chart that shows the top 10 most frew codons for influenza
#c. compare the two results and show the most freq codons aove the two genome
#d. show in the output of the console the top 3 aminoacids of each genome
#e. formulate a prompt for the AI such that the 3 aminoacids are used to ask the AI which 
#food contain those aminoacids 


import os
from collections import Counter
import matplotlib.pyplot as plt


def read_fasta(path: str) -> str:
    """Read first sequence from FASTA file."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing file: {path}")
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    return "".join(seq).upper()


def count_codons(seq: str) -> Counter:
    """Count codon frequencies (frame 0)."""
    c = Counter()
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if "N" in codon or len(codon) < 3:
            continue
        c[codon] += 1
    return c


GENETIC_CODE = {
    "UUU":"Phe","UUC":"Phe","UUA":"Leu","UUG":"Leu",
    "UCU":"Ser","UCC":"Ser","UCA":"Ser","UCG":"Ser",
    "UAU":"Tyr","UAC":"Tyr","UAA":"Stop","UAG":"Stop",
    "UGU":"Cys","UGC":"Cys","UGA":"Stop","UGG":"Trp",
    "CUU":"Leu","CUC":"Leu","CUA":"Leu","CUG":"Leu",
    "CCU":"Pro","CCC":"Pro","CCA":"Pro","CCG":"Pro",
    "CAU":"His","CAC":"His","CAA":"Gln","CAG":"Gln",
    "CGU":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg",
    "AUU":"Ile","AUC":"Ile","AUA":"Ile","AUG":"Met",
    "ACU":"Thr","ACC":"Thr","ACA":"Thr","ACG":"Thr",
    "AAU":"Asn","AAC":"Asn","AAA":"Lys","AAG":"Lys",
    "AGU":"Ser","AGC":"Ser","AGA":"Arg","AGG":"Arg",
    "GUU":"Val","GUC":"Val","GUA":"Val","GUG":"Val",
    "GCU":"Ala","GCC":"Ala","GCA":"Ala","GCG":"Ala",
    "GAU":"Asp","GAC":"Asp","GAA":"Glu","GAG":"Glu",
    "GGU":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly",
}

def codon_to_amino(codon: str) -> str:
    """Convert DNA codon (T->U) to amino acid name."""
    return GENETIC_CODE.get(codon.replace("T","U"), "?")


def top3_amino(counter: Counter):
    aa_counts = Counter()
    for codon, n in counter.items():
        aa = codon_to_amino(codon)
        if aa not in ("Stop", "?"):
            aa_counts[aa] += n
    return aa_counts.most_common(3)



if __name__ == "__main__":
    covid_seq = read_fasta("covid.fasta")
    influenza_seq = read_fasta("influenza.fasta")

    covid_codons = count_codons(covid_seq)
    influenza_codons = count_codons(influenza_seq)

    # a+b)
    covid_top10 = covid_codons.most_common(10)
    influenza_top10 = influenza_codons.most_common(10)

    covid_labels, covid_values = zip(*covid_top10)
    flu_labels, flu_values = zip(*influenza_top10)

    x = range(len(covid_labels))
    plt.figure(figsize=(10,6))
    plt.bar(x, covid_values, width=0.4, label="COVID-19", align="center")
    plt.bar([i + 0.4 for i in x], flu_values, width=0.4, label="Influenza", align="center")
    plt.xticks([i + 0.2 for i in x], covid_labels, rotation=45)
    plt.title("Top 10 Most Frequent Codons — COVID-19 vs Influenza")
    plt.xlabel("Codon")
    plt.ylabel("Frequency")
    plt.legend()
    plt.tight_layout()
    plt.show()

    # c) 
    covid_max = max(covid_codons.items(), key=lambda x: x[1])
    flu_max = max(influenza_codons.items(), key=lambda x: x[1])

    if covid_max[1] > flu_max[1]:
        most_frequent = ("COVID-19", covid_max[0], covid_max[1])
    else:
        most_frequent = ("Influenza", flu_max[0], flu_max[1])

    print(f"\nMost frequent codon overall: {most_frequent[1]} "
          f"({most_frequent[2]} occurrences) in {most_frequent[0]} genome.")

    # d) 
    covid_top3 = top3_amino(covid_codons)
    flu_top3 = top3_amino(influenza_codons)

    print("\nTop 3 amino acids in COVID-19:")
    for aa, count in covid_top3:
        print(f"  {aa}: {count}")

    print("\nTop 3 amino acids in Influenza:")
    for aa, count in flu_top3:
        print(f"  {aa}: {count}")



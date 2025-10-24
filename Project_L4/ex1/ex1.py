#Implement an application that converts the coding region of a gene into an amino acid sequence.
# Use the genetic code from from below. Take a seq od DNA and convert it into ARN seq.

CODON_TABLE = {
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

def dna_to_rna(dna: str) -> str:
    """Convert DNA sequence into RNA sequence (T -> U)."""
    dna = dna.upper().replace(" ", "")
    valid = "".join(base for base in dna if base in "ACGT")
    return valid.replace("T", "U")

def translate_rna(rna: str) -> str:
    """Translate RNA sequence into amino acid sequence."""
    if not rna:
        return ""
    aa_seq = []
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i+3]
        if len(codon) < 3:
            break
        amino = CODON_TABLE.get(codon, "")
        if amino == "Stop" or amino == "":
            break
        aa_seq.append(amino)
    return "-".join(aa_seq)


if __name__ == "__main__":
    dna_seq = input("Enter DNA sequence: ").strip()
    rna_seq = dna_to_rna(dna_seq)
    protein_seq = translate_rna(rna_seq)
    print("RNA sequence:", rna_seq)
    print("Amino acid sequence:", protein_seq)


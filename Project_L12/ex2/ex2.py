import os
import matplotlib.pyplot as plt

MOTIF_TO_FIND = "AATAAA" 
FOLDER_NAME = "influenza"

def read_genome(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    sequence = ""
    for line in lines:
        if not line.startswith(">"):
            sequence += line.strip().upper()
    return sequence

def calculate_signal(genome, motif):
    motif_length = len(motif)
    signal = []
    
    for i in range(len(genome) - motif_length + 1):
        window_sequence = genome[i : i + motif_length]
        
        score = 0
        for j in range(motif_length):
            if window_sequence[j] == motif[j]:
                score += 1
        
        signal.append(score)
        
    return signal

def generate_chart(filename, signal, motif_length):
    plt.figure(figsize=(10, 4))
    plt.plot(signal, color='blue', linewidth=0.8, label='Signal Strength')
    
    plt.axhline(y=motif_length, color='red', linestyle='--', label='Perfect Match')
    
    plt.title(f"Motif Signal '{MOTIF_TO_FIND}' in: {filename}")
    plt.xlabel("Position in Genome")
    plt.ylabel("Match Score")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    print(f"    [Displaying chart for {filename} - Close window to continue...]")
    plt.show()

def main():
    current_path = os.getcwd()
    full_folder_path = os.path.join(current_path, FOLDER_NAME)

    if not os.path.exists(full_folder_path):
        print(f"ERROR: Could not find folder '{FOLDER_NAME}' at path: {full_folder_path}")
        return

    files = [f for f in os.listdir(full_folder_path) if f.endswith(('.fasta', '.fna', '.txt'))]
    
    if not files:
        print(f"No genome files found in '{FOLDER_NAME}'.")
        return

    print(f"Found {len(files)} genome files. Starting analysis...\n")

    for filename in files:
        file_path = os.path.join(full_folder_path, filename)
        
        print(f"Processing: {filename}...")
        
        sequence = read_genome(file_path)
        
        signal = calculate_signal(sequence, MOTIF_TO_FIND)
        
        print(f"    Genome Length: {len(sequence)} bp")
        generate_chart(filename, signal, len(MOTIF_TO_FIND))


if __name__ == "__main__":
    main()
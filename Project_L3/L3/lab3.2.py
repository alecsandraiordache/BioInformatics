#Use the AI to design an app that uses the sliding window methodology in order to scan a DNA seq 
#from a FASTA FILE, and display the melting temp along the seq, by using a chart. the chart must
#have 2 signals, one for each formula. 
#note: the sliding window ahould have 9 positions

import math
import numpy as np
import matplotlib.pyplot as plt
from tkinter import Tk, Button, filedialog, messagebox

WINDOW_SIZE = 9
NA_MOLAR = 0.05 


def read_first_fasta_sequence(path: str) -> str:
    seq_parts = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_parts:
                    break
                continue
            seq_parts.append(line)
    if not seq_parts:
        raise ValueError("No sequence found in FASTA file.")
    return "".join(seq_parts).upper().replace(" ", "").replace("\t", "")


def tm_simple(window: str) -> float:
    A = window.count("A")
    T = window.count("T")
    G = window.count("G")
    C = window.count("C")
    return 4 * (G + C) + 2 * (A + T)


def tm_salt_adjusted(window: str, na_molar: float = NA_MOLAR) -> float:
    L = len(window)
    gc_percent = 100.0 * (window.count("G") + window.count("C")) / L
    return 81.5 + 16.6 * math.log10(na_molar) + 0.41 * gc_percent - 600.0 / L


def smooth_interpolate(x, y, num_points=800):
    """Smooth curve with cubic spline if available."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_dense = np.linspace(x.min(), x.max(), num_points)
    try:
        from scipy.interpolate import make_interp_spline
        y_dense = make_interp_spline(x, y, k=3)(x_dense)
    except Exception:
        y_dense = np.interp(x_dense, x, y)
    return x_dense, y_dense



def process_fasta_file(path: str):
    seq = read_first_fasta_sequence(path)
    if len(seq) < WINDOW_SIZE:
        messagebox.showerror("Error", f"Sequence shorter than window size ({WINDOW_SIZE})")
        return

    x = []
    y_simple = []
    y_salt = []

    for i in range(len(seq) - WINDOW_SIZE + 1):
        w = seq[i:i + WINDOW_SIZE]
        x.append(i + 1)
        y_simple.append(tm_simple(w))
        y_salt.append(tm_salt_adjusted(w))

    x_s, y_simple_s = smooth_interpolate(x, y_simple)
    _,   y_salt_s   = smooth_interpolate(x, y_salt)

    plt.figure()
    plt.plot(x_s, y_simple_s, label="Tm (simple formula)", linewidth=2)
    plt.plot(x_s, y_salt_s,   label="Tm (salt-adjusted)", linewidth=2)
    plt.xlabel(f"Window start position (window size = {WINDOW_SIZE})")
    plt.ylabel("Melting temperature (°C)")
    plt.title("Melting Temperature Along DNA Sequence")
    plt.legend()
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.tight_layout()
    plt.show()


def open_fasta():
    path = filedialog.askopenfilename(
        title="Select FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna *.txt"), ("All files", "*.*")]
    )
    if path:
        try:
            process_fasta_file(path)
        except Exception as e:
            messagebox.showerror("Error", str(e))


def main():
    root = Tk()
    root.title("DNA Melting Temperature (Sliding Window)")
    root.geometry("400x150")

    btn = Button(root, text="Open FASTA File", command=open_fasta,
                 font=("Arial", 14), bg="#4CAF50", fg="white", width=20)
    btn.pack(pady=40)

    root.mainloop()


if __name__ == "__main__":
    main()

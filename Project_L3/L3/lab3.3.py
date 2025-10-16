#Show the min and the maax values of the 2 signals. Also, allow the user to set
# the trashold (like a filter) that is able to take it in consideration only the values above the 
#trashold. these values should be presented to the user in a second chart as a horizontal bars
#thus the chunks of the signal that are above the 
#signal are shown as a horizontal bar over the signals. Wherever the signal is below the treshole
#that chart should show empty space


import math
import numpy as np
import matplotlib.pyplot as plt
from tkinter import Tk, Button, filedialog, messagebox, simpledialog

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



def plot_signals(x, y1, y2, thr):
    """Plot main chart and threshold-highlight chart."""

    # Smooth signals
    x_s, y1_s = smooth_interpolate(x, y1)
    _,   y2_s = smooth_interpolate(x, y2)


    plt.figure(figsize=(10, 7))
    plt.subplot(2, 1, 1)
    plt.plot(x_s, y1_s, label="Tm (simple)", linewidth=2)
    plt.plot(x_s, y2_s, label="Tm (salt-adjusted)", linewidth=2)
    plt.axhline(thr, color="red", linestyle="--", linewidth=1.2, label=f"Threshold = {thr}")
    plt.xlabel(f"Window start (window size = {WINDOW_SIZE})")
    plt.ylabel("Melting Temperature (°C)")
    plt.title("Melting Temperature Along DNA Sequence")
    plt.legend()
    plt.grid(True, linestyle="--", linewidth=0.5)

    
    plt.subplot(2, 1, 2)
    above1 = np.array(y1) > thr
    above2 = np.array(y2) > thr
    plt.bar(x, above1.astype(int)*1.0, width=1.0, color="blue", alpha=0.6, label="Simple > threshold")
    plt.bar(x, above2.astype(int)*-1.0, width=1.0, color="orange", alpha=0.6, label="Salt > threshold")

    plt.axhline(0, color="black", linewidth=1)
    plt.yticks([-1, 0, 1], ["Salt > thr", "", "Simple > thr"])
    plt.xlabel("Window start")
    plt.title("Regions Above Threshold")
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.show()




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


    msg = (
        f"Tm (simple):  min={min(y_simple):.2f} °C   max={max(y_simple):.2f} °C\n"
        f"Tm (salt-adjusted):  min={min(y_salt):.2f} °C   max={max(y_salt):.2f} °C"
    )
    messagebox.showinfo("Signal Information", msg)

    thr = simpledialog.askfloat("Threshold", "Enter threshold value (°C):", minvalue=0.0, initialvalue=50.0)
    if thr is None:
        return

    plot_signals(x, y_simple, y_salt, thr)




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

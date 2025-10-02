#https://ftp.ncbi.nlm.nih.gov/genomes/HUMAN_MICROBIOM/Bacteria/Bacteroides_3_1_33FAA_uid38353/
import re
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from tkinter.scrolledtext import ScrolledText


# =============================
# Data structures & FASTA utils
# =============================
@dataclass
class FastaRecord:
    header: str
    seq: str
    seq_id: str
    species: Optional[str]
    description: str


def parse_fasta(path: Path) -> List[FastaRecord]:
    """Parse a FASTA file into records.

    - Expects first line of each record to start with '>'.
    - Allows sequences on subsequent lines (multiple lines) until next header.
    - Trims whitespace; ignores empty lines between sequences.
    """
    records: List[FastaRecord] = []

    with path.open("r", encoding="utf-8") as f:
        header = None
        seq_lines = []
        for line_num, raw in enumerate(f, start=1):
            line = raw.rstrip("\n")
            if not line:
                # Skip stray blank lines
                continue
            if line.startswith(">"):
                # Commit previous record
                if header is not None:
                    records.append(_build_record(header, seq_lines))
                header = line[1:].strip()
                seq_lines = []
            else:
                if header is None:
                    raise ValueError(
                        f"FASTA format error at line {line_num}: sequence before first header.")
                seq_lines.append(line.strip())
        # Last record
        if header is not None:
            records.append(_build_record(header, seq_lines))

    if not records:
        raise ValueError("No FASTA records found.")

    return records


def _build_record(header: str, seq_lines: List[str]) -> FastaRecord:
    seq = "".join(seq_lines).replace(" ", "")
    seq_id, species, desc = _parse_header(header)
    return FastaRecord(header=header, seq=seq, seq_id=seq_id, species=species, description=desc)


def _parse_header(header: str):
    """Heuristically parse ID, species, description from a FASTA header.

    Rules supported:
    - First token (space-delimited) is treated as ID (common practice)
    - Species captured from common annotations like [organism=...], OS=..., or sp=...
    - If absent, try to guess a binomial like 'Homo sapiens' or 'Escherichia coli'
    """
    tokens = header.split()
    seq_id = tokens[0] if tokens else "seq1"

    # Try common patterns
    species = None
    m = re.search(r"\[organism=([^\]]+)\]", header, flags=re.IGNORECASE)
    if m:
        species = m.group(1).strip()
    if species is None:
        m = re.search(r"(?:OS|organism|sp|species)=([^|]+?)(?:\||\s|$)", header, flags=re.IGNORECASE)
        if m:
            species = m.group(1).strip()

    # Fallback: look for two capitalized words (very heuristic)
    if species is None:
        m = re.search(r"\b([A-Z][a-z]+\s+[a-z]+)\b", header)
        if m:
            species = m.group(1)

    # Description is everything after the first token
    desc = header[len(seq_id):].strip() if len(header) > len(seq_id) else ""

    return seq_id, species, desc


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(textwrap.wrap(seq, width=width))


def detect_seq_type(seq: str) -> str:
    """Naively detect sequence type: DNA, RNA, or Protein."""
    s = seq.upper()
    # If U present and T absent, likely RNA
    if "U" in s and "T" not in s and re.fullmatch(r"[ACGUNX\-]+", s) is not None:
        return "RNA"
    # If only DNA alphabet
    if re.fullmatch(r"[ACGTNX\-]+", s) is not None:
        return "DNA"
    # Otherwise protein (very rough)
    return "Protein"


# ===========================================
# Assignment 1 & 2 integration (placeholders)
# ===========================================
# Replace the bodies of these functions with the logic from your Assignment 1 and 2.
# The GUI will pass the selected sequence string "seq" and display your returned result.

def assignment1_algorithm(seq: str) -> str:
    """Demo implementation: basic statistics (GC content for nucleic acids; length otherwise)."""
    seq_type = detect_seq_type(seq)
    n = len(seq)
    if seq_type in ("DNA", "RNA"):
        s = seq.upper()
        g = s.count("G")
        c = s.count("C")
        gc = (g + c) / n * 100 if n else 0.0
        return f"Assignment 1 (demo): length={n} bp, GC%={gc:.2f}"
    else:
        return f"Assignment 1 (demo): length={n} aa"


def assignment2_algorithm(seq: str) -> str:
    """Demo implementation: simple composition table (top 5 symbols)."""
    from collections import Counter
    c = Counter(seq.upper())
    items = sorted(c.items(), key=lambda kv: (-kv[1], kv[0]))
    top = ", ".join(f"{k}:{v}" for k, v in items[:5])
    return f"Assignment 2 (demo): most frequent symbols → {top}"


# =====================
# GUI implementation
# =====================
class FastaGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FASTA Analyzer – Assignment 1 & 2 Integrator")
        self.geometry("980x700")
        self.minsize(900, 640)

        self.records: List[FastaRecord] = []
        self.current_idx: Optional[int] = None

        self._build_widgets()

    def _build_widgets(self):
        # Top frame: file controls
        top = ttk.Frame(self)
        top.pack(fill=tk.X, padx=12, pady=8)

        self.btn_open = ttk.Button(top, text="Choose FASTA file…", command=self.open_fasta)
        self.btn_open.pack(side=tk.LEFT)

        self.file_label_var = tk.StringVar(value="No file loaded")
        ttk.Label(top, textvariable=self.file_label_var).pack(side=tk.LEFT, padx=10)

        # Middle frame: record selector + metadata
        mid = ttk.Frame(self)
        mid.pack(fill=tk.X, padx=12, pady=4)

        ttk.Label(mid, text="Sequence:").grid(row=0, column=0, sticky=tk.W, padx=(0, 6))
        self.combo_records = ttk.Combobox(mid, state="readonly", width=50)
        self.combo_records.bind("<<ComboboxSelected>>", self._on_record_change)
        self.combo_records.grid(row=0, column=1, sticky=tk.W)

        # Metadata grid
        self.id_var = tk.StringVar()
        self.species_var = tk.StringVar()
        self.type_var = tk.StringVar()
        self.length_var = tk.StringVar()

        row2 = 1
        ttk.Label(mid, text="ID:").grid(row=row2, column=0, sticky=tk.E, padx=(0, 6), pady=2)
        ttk.Label(mid, textvariable=self.id_var).grid(row=row2, column=1, sticky=tk.W, pady=2)

        row2 += 1
        ttk.Label(mid, text="Species:").grid(row=row2, column=0, sticky=tk.E, padx=(0, 6), pady=2)
        ttk.Label(mid, textvariable=self.species_var).grid(row=row2, column=1, sticky=tk.W, pady=2)

        row2 += 1
        ttk.Label(mid, text="Type:").grid(row=row2, column=0, sticky=tk.E, padx=(0, 6), pady=2)
        ttk.Label(mid, textvariable=self.type_var).grid(row=row2, column=1, sticky=tk.W, pady=2)

        row2 += 1
        ttk.Label(mid, text="Length:").grid(row=row2, column=0, sticky=tk.E, padx=(0, 6), pady=2)
        ttk.Label(mid, textvariable=self.length_var).grid(row=row2, column=1, sticky=tk.W, pady=2)

        # Sequence view & actions
        body = ttk.Panedwindow(self, orient=tk.HORIZONTAL)
        body.pack(fill=tk.BOTH, expand=True, padx=12, pady=8)

        left = ttk.Frame(body)
        right = ttk.Frame(body)
        body.add(left, weight=3)
        body.add(right, weight=2)

        # Left pane: wrapped sequence
        ttk.Label(left, text="Sequence (wrapped to 80 characters):").pack(anchor=tk.W)
        self.seq_text = ScrolledText(left, height=20, wrap=tk.NONE, font=("Courier New", 10))
        self.seq_text.pack(fill=tk.BOTH, expand=True)
        self.seq_text.configure(state=tk.DISABLED)

        action_bar = ttk.Frame(left)
        action_bar.pack(fill=tk.X, pady=(6, 0))
        ttk.Button(action_bar, text="Copy sequence", command=self.copy_sequence).pack(side=tk.LEFT)
        ttk.Button(action_bar, text="Save wrapped as FASTA…", command=self.save_wrapped_fasta).pack(side=tk.LEFT, padx=8)

        # Right pane: run assignments and see results
        ttk.Label(right, text="Run assignments:").pack(anchor=tk.W)
        btns = ttk.Frame(right)
        btns.pack(fill=tk.X, pady=4)
        ttk.Button(btns, text="Run Assignment 1", command=self.run_a1).pack(side=tk.LEFT)
        ttk.Button(btns, text="Run Assignment 2", command=self.run_a2).pack(side=tk.LEFT, padx=8)

        ttk.Label(right, text="Output:").pack(anchor=tk.W, pady=(8, 0))
        self.output = ScrolledText(right, height=18, wrap=tk.WORD, font=("Segoe UI", 10))
        self.output.pack(fill=tk.BOTH, expand=True)
        self.output.configure(state=tk.DISABLED)

        # Status bar
        self.status_var = tk.StringVar(value="Welcome! Load a FASTA file to begin.")
        status = ttk.Label(self, textvariable=self.status_var, relief=tk.SUNKEN, anchor=tk.W)
        status.pack(fill=tk.X, side=tk.BOTTOM)

    # ---------------
    # Event handlers
    # ---------------
    def open_fasta(self):
        path_str = filedialog.askopenfilename(
            title="Open FASTA",
            filetypes=[
                ("FASTA files", "*.fa *.fasta *.fna *.ffn *.faa *.frn"),
                ("All files", "*.*"),
            ],
        )
        if not path_str:
            return

        path = Path(path_str)
        try:
            self.records = parse_fasta(path)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to parse FASTA file:\n{e}")
            return

        self.file_label_var.set(str(path.name))
        self.status_var.set(f"Loaded {len(self.records)} record(s) from {path.name}")

        # Populate selector
        self.combo_records["values"] = [self._record_label(i, r) for i, r in enumerate(self.records, start=1)]
        if self.records:
            self.combo_records.current(0)
            self._load_record(0)

    def _record_label(self, idx: int, rec: FastaRecord) -> str:
        sp = f" | {rec.species}" if rec.species else ""
        return f"{idx}. {rec.seq_id}{sp} ({len(rec.seq)})"

    def _on_record_change(self, _event=None):
        i = self.combo_records.current()
        if i >= 0:
            self._load_record(i)

    def _load_record(self, idx: int):
        self.current_idx = idx
        rec = self.records[idx]
        self.id_var.set(rec.seq_id)
        self.species_var.set(rec.species or "—")
        self.type_var.set(detect_seq_type(rec.seq))
        self.length_var.set(f"{len(rec.seq)}")

        wrapped = wrap_sequence(rec.seq, 80)
        self.seq_text.configure(state=tk.NORMAL)
        self.seq_text.delete("1.0", tk.END)
        self.seq_text.insert("1.0", wrapped)
        self.seq_text.configure(state=tk.DISABLED)

        self._write_output(f"Loaded sequence: {rec.seq_id}\nHeader: >{rec.header}\n")

    def _write_output(self, text: str, clear: bool = False):
        self.output.configure(state=tk.NORMAL)
        if clear:
            self.output.delete("1.0", tk.END)
        self.output.insert(tk.END, text + ("\n" if not text.endswith("\n") else ""))
        self.output.see(tk.END)
        self.output.configure(state=tk.DISABLED)

    def copy_sequence(self):
        if self.current_idx is None:
            return
        rec = self.records[self.current_idx]
        self.clipboard_clear()
        self.clipboard_append(rec.seq)
        self.status_var.set("Sequence copied to clipboard.")

    def save_wrapped_fasta(self):
        if self.current_idx is None:
            return
        rec = self.records[self.current_idx]
        save_path = filedialog.asksaveasfilename(
            title="Save FASTA",
            defaultextension=".fasta",
            initialfile=f"{rec.seq_id}.fasta",
            filetypes=[("FASTA", "*.fasta *.fa"), ("All files", "*.*")],
        )
        if not save_path:
            return

        with open(save_path, "w", encoding="utf-8") as out:
            header = rec.header or rec.seq_id
            out.write(">" + header + "\n")
            out.write(wrap_sequence(rec.seq, 80) + "\n")
        self.status_var.set(f"Saved FASTA to {Path(save_path).name}")

    def _get_current_seq(self) -> Optional[str]:
        if self.current_idx is None:
            messagebox.showinfo("No sequence", "Please load a FASTA file and select a sequence.")
            return None
        return self.records[self.current_idx].seq

    def run_a1(self):
        seq = self._get_current_seq()
        if not seq:
            return
        try:
            result = assignment1_algorithm(seq)
        except Exception as e:
            messagebox.showerror("Assignment 1 error", str(e))
            return
        self._write_output(f"[Assignment 1] {result}")

    def run_a2(self):
        seq = self._get_current_seq()
        if not seq:
            return
        try:
            result = assignment2_algorithm(seq)
        except Exception as e:
            messagebox.showerror("Assignment 2 error", str(e))
            return
        self._write_output(f"[Assignment 2] {result}")


# --------------
# Entry point
# --------------
if __name__ == "__main__":
    try:
        app = FastaGUI()
        app.mainloop()
    except Exception as exc:
        messagebox.showerror("Fatal error", str(exc))
        sys.exit(1)

import matplotlib.pyplot as plt
import statistics

S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"

def cg_percent(seq: str) -> float:
    seq = seq.upper()
    L = len(seq)
    if L == 0:
        return 0.0
    cg = sum(1 for b in seq if b in ("C", "G"))
    return round((cg / L) * 100.0, 2)


def _kappa_raw(seq: str) -> float:
    A = seq.upper()
    L = len(A)
    if L < 2:
        return 0.0

    N = L - 1
    T = 0.0

    for u in range(1, N + 1):
        matches = 0
        B_len = L - u
        for i in range(B_len):
            if A[i] == A[u + i]:
                matches += 1
        T += (matches / B_len) * 100.0

    return T / N


_KAPPA_SCALE = 27.53 / _kappa_raw(S)


def kappa_ic(seq: str) -> float:
    raw = _kappa_raw(seq)
    return round(raw * _KAPPA_SCALE, 2)


def sliding_windows(seq: str, win_size: int):
    seq = seq.upper()
    return [seq[i:i+win_size] for i in range(len(seq) - win_size + 1)]


def compute_pattern(seq: str, win_size: int):
    windows = sliding_windows(seq, win_size)
    xs = [cg_percent(w) for w in windows]
    ys = [kappa_ic(w) for w in windows]
    return xs, ys


def center_of_weight(xs, ys):
    if not xs or not ys:
        return 0.0, 0.0
    return statistics.mean(xs), statistics.mean(ys)


def main():
    CG = cg_percent(S)
    IC = kappa_ic(S)

    print(f"CG% for S = {CG:.2f}")     
    print(f"Kappa IC for S = {IC:.2f}") 

    win_size = 30

    xs, ys = compute_pattern(S, win_size)

    cx, cy = center_of_weight(xs, ys)
    print(f"Center of weight of pattern: ({cx:.2f}, {cy:.2f})")

    plt.figure()
    plt.scatter(xs, ys, s=20)
    plt.xlabel("C+G%")
    plt.ylabel("Kappa IC")
    plt.title("DNA Pattern (Sliding windows of 30 bp)")
    plt.grid(True)

    plt.figure()
    plt.scatter([cx], [cy], color="red", s=80)
    plt.xlabel("C+G% (center)")
    plt.ylabel("IC (center)")
    plt.title("Center of Weight of the Pattern")
    plt.grid(True)

    plt.show()


if __name__ == "__main__":
    main()

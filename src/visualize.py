"""
Visualize the Needleman-Wunsch DP matrix as a heatmap.
Requires: matplotlib, numpy
Run: python src/visualize.py
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
except ImportError:
    print("Install matplotlib: pip install matplotlib numpy")
    sys.exit(1)

from src.nw import needleman_wunsch


def plot_dp_matrix(seq1: str, seq2: str, match=1, mismatch=-1, gap=-1):
    a1, a2, score, dp = needleman_wunsch(seq1, seq2, match, mismatch, gap)
    mat = np.array(dp)

    fig, ax = plt.subplots(figsize=(max(6, len(seq2) * 0.7 + 1),
                                     max(5, len(seq1) * 0.7 + 1)))

    im = ax.imshow(mat, cmap="RdYlGn", aspect="auto")
    plt.colorbar(im, ax=ax, label="Score")

    # Labels
    ax.set_xticks(range(len(seq2) + 1))
    ax.set_xticklabels(["-"] + list(seq2), fontsize=12, fontfamily="monospace")
    ax.set_yticks(range(len(seq1) + 1))
    ax.set_yticklabels(["-"] + list(seq1), fontsize=12, fontfamily="monospace")

    # Cell values
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, str(mat[i, j]), ha="center", va="center",
                    fontsize=9, color="black", fontfamily="monospace")

    ax.set_title(
        f"NW DP Matrix  |  score={score}\n"
        f"{seq1}  ↔  {seq2}  "
        f"[match={match} mismatch={mismatch} gap={gap}]",
        fontsize=11, pad=12
    )
    ax.xaxis.set_label_position("top")
    ax.xaxis.tick_top()

    plt.tight_layout()

    out = Path(__file__).parent.parent / "dp_matrix.png"
    plt.savefig(out, dpi=150)
    print(f"  Saved → {out}")
    plt.show()


if __name__ == "__main__":
    seq1 = sys.argv[1] if len(sys.argv) > 2 else "GATTACA"
    seq2 = sys.argv[2] if len(sys.argv) > 2 else "GCATGCU"
    plot_dp_matrix(seq1, seq2)

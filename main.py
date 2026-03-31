#!/usr/bin/env python3
"""
Needleman-Wunsch Global Sequence Aligner
Usage:
    python main.py                          # interactive
    python main.py seq1.txt seq2.txt        # from files
    python main.py --match 2 --gap -2       # custom scoring
"""

import sys
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from src.nw import needleman_wunsch, print_matrix


def load_sequence(path_or_seq: str) -> str:
    """Accept raw sequence string or path to .txt / .fasta file."""
    p = Path(path_or_seq)
    if p.exists():
        text = p.read_text().strip()
        if text.startswith(">"):                # FASTA
            lines = text.splitlines()
            return "".join(l for l in lines if not l.startswith(">")).upper()
        return text.upper()
    return path_or_seq.upper()


def format_alignment(a1: str, a2: str) -> str:
    mid = "".join("|" if c1 == c2 and c1 != "-" else " " for c1, c2 in zip(a1, a2))
    return f"  {a1}\n  {mid}\n  {a2}"


def run(seq1: str, seq2: str, match: int, mismatch: int, gap: int, show_matrix: bool):
    print(f"\n  seq1 ({len(seq1)} bp): {seq1}")
    print(f"  seq2 ({len(seq2)} bp): {seq2}")
    print(f"  scoring  match={match}  mismatch={mismatch}  gap={gap}\n")

    a1, a2, score, dp = needleman_wunsch(seq1, seq2, match, mismatch, gap)

    if show_matrix:
        print("  ── DP matrix ──")
        print_matrix(dp, seq1, seq2)
        print()

    print("  ── Alignment ──")
    print(format_alignment(a1, a2))
    print(f"\n  Score: {score}")

    # Optional: compare with Biopython if installed
    try:
        from Bio import pairwise2
        bio = pairwise2.align.globalms(seq1, seq2, match, mismatch, gap, gap)
        bio_score = bio[0].score
        match_str = "✓ match" if bio_score == score else f"✗ mismatch (Bio={bio_score})"
        print(f"  Biopython check: {match_str}")
    except ImportError:
        pass


def main():
    parser = argparse.ArgumentParser(description="Needleman-Wunsch aligner")
    parser.add_argument("seq1", nargs="?", help="sequence or .txt/.fasta path")
    parser.add_argument("seq2", nargs="?", help="sequence or .txt/.fasta path")
    parser.add_argument("--match",    type=int, default=1,  help="match score (default 1)")
    parser.add_argument("--mismatch", type=int, default=-1, help="mismatch penalty (default -1)")
    parser.add_argument("--gap",      type=int, default=-1, help="gap penalty (default -1)")
    parser.add_argument("--matrix",   action="store_true",  help="print DP matrix")
    args = parser.parse_args()

    if args.seq1 and args.seq2:
        s1 = load_sequence(args.seq1)
        s2 = load_sequence(args.seq2)
    else:
        print("  Needleman-Wunsch Aligner — interactive mode")
        s1 = input("  seq1: ").strip().upper()
        s2 = input("  seq2: ").strip().upper()

    run(s1, s2, args.match, args.mismatch, args.gap, args.matrix)


if __name__ == "__main__":
    main()

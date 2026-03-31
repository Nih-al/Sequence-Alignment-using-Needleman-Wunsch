"""
Tests for Needleman-Wunsch implementation.
Run: python -m pytest tests/ -v
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.nw import needleman_wunsch


def test_classic():
    a1, a2, score, _ = needleman_wunsch("GATTACA", "GCATGCU")
    assert score == 0, f"Expected 0, got {score}"

def test_identical():
    a1, a2, score, _ = needleman_wunsch("ACGT", "ACGT")
    assert a1 == a2 == "ACGT"
    assert score == 4

def test_all_gaps():
    a1, a2, score, _ = needleman_wunsch("AAAA", "CCCC")
    assert score == -4   # 4 mismatches at -1 each

def test_empty_seq1():
    a1, a2, score, _ = needleman_wunsch("", "ACG")
    assert a1 == "---"
    assert a2 == "ACG"
    assert score == -3

def test_empty_seq2():
    a1, a2, score, _ = needleman_wunsch("ACG", "")
    assert a2 == "---"
    assert a1 == "ACG"
    assert score == -3

def test_single_char_match():
    a1, a2, score, _ = needleman_wunsch("A", "A")
    assert score == 1

def test_single_char_mismatch():
    a1, a2, score, _ = needleman_wunsch("A", "C")
    assert score == -1

def test_custom_scoring():
    _, _, score, _ = needleman_wunsch("ACGT", "ACGT", match=2, mismatch=-2, gap=-2)
    assert score == 8

def test_alignment_length():
    """Aligned sequences must always have the same length."""
    a1, a2, _, _ = needleman_wunsch("AGCT", "TAGCTA")
    assert len(a1) == len(a2)

def test_no_stray_bases():
    """All original bases must appear in alignment (gaps notwithstanding)."""
    seq1, seq2 = "AGTACG", "ACATCG"
    a1, a2, _, _ = needleman_wunsch(seq1, seq2)
    assert a1.replace("-", "") == seq1
    assert a2.replace("-", "") == seq2


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = failed = 0
    for t in tests:
        try:
            t()
            print(f"  PASS  {t.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL  {t.__name__}  →  {e}")
            failed += 1
    print(f"\n  {passed}/{passed+failed} passed")

# Needleman-Wunsch Global Sequence Aligner

Manual implementation of the NW algorithm with DP matrix visualization and Biopython validation.

## Structure

```
nw_project/
├── main.py              # CLI entrypoint
├── src/
│   ├── nw.py            # Core algorithm (DP + traceback)
│   └── visualize.py     # Heatmap visualization
├── tests/
│   └── test_nw.py       # 10 edge-case tests
└── data/
    ├── seq1.fasta
    └── seq2.fasta
```

## Usage

```bash
# Interactive
python main.py

# Positional args (raw strings)
python main.py GATTACA GCATGCU

# From FASTA files
python main.py data/seq1.fasta data/seq2.fasta

# Custom scoring + show DP matrix
python main.py GATTACA GCATGCU --match 2 --mismatch -2 --gap -2 --matrix
```

## Tests

```bash
python tests/test_nw.py          # direct
python -m pytest tests/ -v       # with pytest
```

## Visualization

```bash
python src/visualize.py GATTACA GCATGCU
# → saves dp_matrix.png
```

## Algorithm

Dynamic programming with O(nm) time and space.

Recurrence:
```
dp[i][j] = max(
    dp[i-1][j-1] + match/mismatch,
    dp[i-1][j]   + gap,
    dp[i][j-1]   + gap
)
```

Traceback reconstructs alignment by walking backwards from `dp[n][m]`.

## Dependencies

```
biopython   # optional, for score validation
matplotlib  # optional, for heatmap
numpy       # optional, for heatmap
```

Install: `pip install biopython matplotlib numpy`

def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    n, m = len(seq1), len(seq2)

    # Initialize DP matrix
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        dp[i][0] = i * gap
    for j in range(m + 1):
        dp[0][j] = j * gap

    # Fill matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = match if seq1[i - 1] == seq2[j - 1] else mismatch
            diag = dp[i - 1][j - 1] + score
            up   = dp[i - 1][j] + gap
            left = dp[i][j - 1] + gap
            dp[i][j] = max(diag, up, left)

    # Traceback
    i, j = n, m
    align1, align2 = "", ""

    while i > 0 and j > 0:
        current = dp[i][j]
        score = match if seq1[i - 1] == seq2[j - 1] else mismatch
        if current == dp[i - 1][j - 1] + score:
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1; j -= 1
        elif current == dp[i - 1][j] + gap:
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    while i > 0:
        align1 = seq1[i - 1] + align1
        align2 = "-" + align2
        i -= 1
    while j > 0:
        align1 = "-" + align1
        align2 = seq2[j - 1] + align2
        j -= 1

    return align1, align2, dp[n][m], dp


def print_matrix(dp, seq1, seq2):
    header = "    -  " + "  ".join(seq2)
    print(header)
    for i, row in enumerate(dp):
        label = "-" if i == 0 else seq1[i - 1]
        print(f"{label}  " + "  ".join(f"{v:2d}" for v in row))

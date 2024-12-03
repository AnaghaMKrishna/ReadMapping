#!/usr/bin/env python3

def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    """
    function to implement NW algorithm

    Args:
        seq_a: first sequence to align
        seq_b: second sequence to align
        match: match score
        mismatch: mismatch penalty
        gap: gap penalty

    Returns:
        A tuple containing aligned sequences and the matched score
    """
    matrix = generate_matrix(len(seq_a), len(seq_b))
    filled_matrix = fill_matrix(matrix = matrix, seq_a = seq_a, seq_b = seq_b, match = match, mismatch = mismatch, gap = gap)
    aln, score = find_optimal_alignment(matrix = filled_matrix, seq_a = seq_a, seq_b = seq_b, match = match, mismatch = mismatch, gap = gap)

    return aln, score

def generate_matrix(row: int, col: int) -> list[list[int]]:
    """
    creates a (row + 1) X (col + 1) matrix. Extra row and column added to fill in initial values.

    Args:
        row: length of seq 1
        col: length of seq 2

    Returns:
        two-dimensional matrix
    """
    mat = [[ 0 for j in range(col + 1)] for i in range(row + 1)]
    return mat

def fill_matrix(matrix: list[list[int]], seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> list[list[int]]:
    """
    fill the matrix with scores based on the sequence alignment

    Args:
        matrix: 2D initialized matrix
        seq_a: first sequence
        seq_b: second sequence
        match: match score
        mismatch: mismatch penalty
        gap: gap penalty

    Returns:
        score matrix
    """
    #fill the first row
    for col in range(1, len(seq_b) + 1):
        matrix[0][col] = matrix[0][col - 1] + gap

    #fill the first column
    for row in range(1, len(seq_a) + 1):
        matrix[row][0] = matrix[row - 1][0] + gap

    #fill the rest of the matrix with maximum of diagonal, right or top cell
    for row in range(1, len(seq_a) + 1):
        seq_a_base = seq_a[row - 1]
        for col in range(1, len(seq_b) + 1):
            seq_b_base = seq_b[col - 1]
            matrix[row][col] = max(
                                    (matrix[row - 1][col] + gap), 
                                    (matrix[row][col - 1] + gap), 
                                    (matrix[row - 1][col -1] + match if seq_a_base == seq_b_base else matrix[row - 1][col -1] + mismatch )
                                    )       
    # for row in matrix:
    #     print(row)

    return matrix

def find_optimal_alignment(matrix: list[list[int]], seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    """
    find optimal alignment for the sequences using backtracking the score matrix

    Args:
        matrix: score matrix
        seq_a: first sequence
        seq_b: second sequence
        match: match score
        mismatch: mismatch penalty
        gap: gap penalty

    Returns:
        retruns one optimal alignment of the two sequences
    """
    score = matrix[len(seq_a)][len(seq_b)]

    aln_a = ""
    aln_b = ""
    row = len(seq_a)
    col = len(seq_b)

    #traverse the matrix by starting at bottom right-most cell and reach to (0,0)
    while row > 0 and col > 0:
        #get the index of max value among diagonal, left and top cells
        check_max_list = [matrix[row-1][col-1], matrix[row][col-1], matrix[row-1][col]]
        max_pos = check_max_list.index(max(check_max_list))
        #for match
        if max_pos == 0 and matrix[row][col] == matrix[row-1][col-1] + match:
            aln_b += seq_b[col-1]
            aln_a += seq_a[row-1]
            row -= 1
            col -= 1
        #for gap on seq_b
        elif max_pos == 1 and matrix[row][col] == matrix[row][col-1] + gap:
            aln_a += "-"
            aln_b += seq_b[col-1]
            col -= 1
        #for gap on seq_a
        elif max_pos == 2 and matrix[row][col] == matrix[row-1][col] + gap:
            aln_b += "-"
            aln_a += seq_a[row-1]
            row -= 1
        #mismatch
        else:
            aln_b += seq_b[col-1]
            aln_a += seq_a[row-1]
            row -= 1
            col -= 1
    #reverse the string before returning
    a = aln_a[::-1]
    b = aln_b[::-1]
    return ((a, b), score)

# if __name__ == "__main__":
#     a,b = needleman_wunsch("ATA", "ATCA", 1, -1, -1)
#     print(a)
#     print(b)
#     a,b = needleman_wunsch("TAGTCAT", "TATCAAT", 1, -1, -1)
#     print(a)
#     print(b)
#     a,b = needleman_wunsch("CTTCTCGTCGGTCTCGTGGTTCGGGAAC", "CTTTCATCCACTTCGTTGCCCGGGAAC", 1, -1, -1)
#     print(a)
#     print(b)
import numpy as np


# Parsing the input file to get the 2 sequences in it and store them in 2 string to be further processed
def parse_input(dir1, dir2):
    with open(dir1, 'r') as file1, open(dir2, 'r') as file2:
        sequence1 = file1.read()
        sequence2 = file2.read()

    return sequence1, sequence2


def initialize(seq1, seq2, gap_penalty):
    width = len(seq1) + 1
    height = len(seq2) + 1

    initialization_matrix = np.zeros((height, width), dtype=int)
    v_arrows = np.zeros((height, width), dtype=bool)
    h_arrows = np.zeros((height, width), dtype=bool)
    d_arrows = np.zeros((height, width), dtype=bool)

    for i in range(height):
        initialization_matrix[i, 0] = gap_penalty * i
        v_arrows[i, 0] = True
    for j in range(width):
        initialization_matrix[0, j] = gap_penalty * j
        h_arrows[0, j] = True

    return initialization_matrix, h_arrows, v_arrows, d_arrows


def matrix_fill(score_matrix, h_arrows, v_arrows, d_arrows, match_score, mismatch_penalty, gap_penalty, seq1, seq2):
    width = score_matrix.shape[0]
    height = score_matrix.shape[1]

    for i in range(1, height):
        for j in range(1, width):
            if seq2[i - 1] == seq1[j - 1]:
                match_flag = True
                mismatch_flag = False
            else:
                mismatch_flag = True
                match_flag = False

            diagonal_score = score_matrix[i - 1, j - 1] + match_score * match_flag + mismatch_penalty * mismatch_flag
            vertical_score = score_matrix[i - 1, j] + gap_penalty
            horizontal_score = score_matrix[i, j - 1] + gap_penalty

            score_matrix[i, j] = max(diagonal_score, vertical_score, horizontal_score)

            if score_matrix[i, j] == diagonal_score:
                d_arrows[i, j] = True
            if score_matrix[i, j] == horizontal_score:
                h_arrows[i, j] = True
            if score_matrix[i, j] == vertical_score:
                v_arrows[i, j] = True

    return score_matrix, h_arrows, v_arrows, d_arrows


def trace_back(score_matrix, h_arrows, v_arrows, d_arrows, seq1, seq2, match_score, mismatch_penalty, gap_penalty):
    width = score_matrix.shape[0]
    height = score_matrix.shape[1]

    i = width - 1
    j = height - 1

    sequence1 = ""
    sequence2 = ""

    while i != 0 and j != 0:

        if d_arrows[i, j]:
            sequence1 += seq1[j-1]
            sequence2 += seq2[i-1]
            i = i - 1
            j = j - 1
        elif h_arrows[i, j]:
            sequence1 += seq1[j-1]
            sequence2 += '-'
            j = j - 1

        elif v_arrows[i, j]:
            sequence1 += '-'
            sequence2 += seq2[i-1]
            i = i - 1

    return sequence1[::-1], sequence2[::-1]  # Return the strings in the correct order


def global_alignment(dir1, dir2, match_score, mismatch_penalty, gap_penalty):
    seq1, seq2 = parse_input(dir1, dir2)
    initialization_matrix, h_arrows, v_arrows, d_arrows = initialize(seq1, seq2, gap_penalty)
    score_matrix, h_arrows, v_arrows, d_arrows = matrix_fill(initialization_matrix, h_arrows, v_arrows, d_arrows,
                                                             match_score, mismatch_penalty, gap_penalty, seq1, seq2)
    possible_layouts = trace_back(score_matrix, h_arrows, v_arrows, d_arrows, seq1, seq2, match_score, mismatch_penalty,
                              gap_penalty)
    # print(score_matrix)
    # print(initialization_matrix)
    # print(h_arrows)
    # print(v_arrows)
    print(possible_layouts)


# Test case for the 2 sequences in lecture 7 with the score from assignment 1
global_alignment('sequence1.txt', 'sequence2.txt', 4, -1, -3)

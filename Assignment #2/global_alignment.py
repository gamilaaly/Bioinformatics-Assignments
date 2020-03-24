import numpy as np


# Parsing the input file to get the 2 sequences in it and store them in 2 string to be further processed
def parse_input(dir1, dir2):
    with open(dir1, 'r') as file1, open(dir2, 'r') as file2:
        sequence1 = file1.read()
        sequence2 = file2.read()

    return sequence1, sequence2


# Returns 2D matrix filled with initial values
# Returns 3D matrix filled with arrows indicates the directions of initial values
def initialize(seq1, seq2, gap_penalty):
    width = len(seq1) + 1
    height = len(seq2) + 1

    # Construct 2D matrix to store the values of best score in each position
    initialization_matrix = np.zeros((height, width), dtype=int)
    # Construct 3D matrix to store the values of direction from which best score was calculated
    arrows_matrix = np.zeros((height, width, 3), dtype=bool)

    # Fill The 1st row and 1st column of the matrix with initial value equals index multiplied by gap penalty
    for i in range(height):
        initialization_matrix[i, 0] = gap_penalty * i
        arrows_matrix[i, 0, 2] = True  # Vertical arrows
    for j in range(width):
        initialization_matrix[0, j] = gap_penalty * j
        arrows_matrix[0, j, 1] = True  # Horizontal arrows

    return initialization_matrix, arrows_matrix


def matrix_fill(score_matrix, arrows_matrix, match_score, mismatch_penalty, gap_penalty, seq1, seq2):
    width = score_matrix.shape[0]
    height = score_matrix.shape[1]

    for j in range(1, height):
        for i in range(1, width):
            # Determine if the 2 sequences have a match in this postion or not, to use it in the maximum score calculation
            if seq2[i - 1] == seq1[j - 1]:
                match_flag = True
                mismatch_flag = False
            else:
                mismatch_flag = True
                match_flag = False

            # Calculation of possible scores for the current position
            diagonal_score = score_matrix[i - 1, j - 1] + match_score * match_flag + mismatch_penalty * mismatch_flag
            vertical_score = score_matrix[i - 1, j] + gap_penalty
            horizontal_score = score_matrix[i, j - 1] + gap_penalty

            # Calculate the best score at each position
            score_matrix[i, j] = max(diagonal_score, vertical_score, horizontal_score)

            # Calculate the direction/s from which the best score is calculated
            if score_matrix[i, j] == diagonal_score:
                arrows_matrix[i, j, 0] = True  # Diagnonal Arrows
            if score_matrix[i, j] == horizontal_score:
                arrows_matrix[i, j, 1] = True  # Horizontal Arrows
            if score_matrix[i, j] == vertical_score:
                arrows_matrix[i, j, 2] = True  # Vertical arrows
    return score_matrix, arrows_matrix


def trace_back(score_matrix, arrows_matrix, seq1, seq2):
    width = score_matrix.shape[0]
    height = score_matrix.shape[1]

    i = width - 1
    j = height - 1

    sequence1 = ""
    sequence2 = ""

    while i != 0 and j != 0:
        if arrows_matrix[i, j, 0]:  # Diagonal Arrow
            sequence1 += seq1[j - 1]
            sequence2 += seq2[i - 1]
            i = i - 1  # Both indices to move along the diagonal
            j = j - 1
        elif arrows_matrix[i, j, 1]:  # Horizontal arrow
            sequence1 += seq1[j - 1]
            sequence2 += '-'  # Add ga
            j = j - 1  # Move Horizontally
        elif arrows_matrix[i, j, 2]:  # Vertical arrow
            sequence1 += '-'
            sequence2 += seq2[i - 1]
            i = i - 1  # Move vertically

    # Invert the sequence to be returned in the correct formula
    sequence1 = sequence1[::-1]
    sequence2 = sequence2[::-1]

    return sequence1, sequence2


# Find all paths to traceback from best score to zero
def find_paths(score_matrix, arrows_matrix, seq1, seq2):
    width = score_matrix.shape[0]
    height = score_matrix.shape[1]

    i = width - 1
    j = height - 1

    paths = []  # Array of indices of positions that have mor than one direction
    possible_layouts = []  # Array of the 2 sequences after alignment

    while i != 0 and j != 0:
        # Determine if there are different baths to take from this position
        if np.count_nonzero(arrows_matrix[i, j, :]) > 1:
            for x in range(3, 0, -1):
                if arrows_matrix[i, j, x - 1]:
                    paths.append([i, j, x - 1]) # Add the different directions that can be taken from this point

        if arrows_matrix[i, j, 0]:  # Diagonal Arrow
            i = i - 1
            j = j - 1
        elif arrows_matrix[i, j, 1]:  # Horizontal arrow
            j = j - 1
        elif arrows_matrix[i, j, 2]:  # Vertical arrow
            i = i - 1

    # If there is only on path, the traceback function should be called only once
    if len(paths) == 0:
        possible_layouts.append(trace_back(score_matrix, arrows_matrix, seq1, seq2))

    # Invert the directions
    paths = paths[::-1]

    # Tracing back All possible paths
    while len(paths) > 0:
        possible_layouts.append(trace_back(score_matrix, arrows_matrix, seq1, seq2))

        if paths[0][2] == 0:  # for the 1st iteration
            arrows_matrix[paths[0][0], paths[0][1], paths[0][2]] = False
            paths.pop(0) # Remove the direction that it's already used
        elif len(paths) == 1: # for the last iteration, There will be only 1 element left in the possible directions to pop out
            arrows_matrix[paths[0][0], paths[0][1], paths[0][2]] = False
            paths.pop(0)
        elif paths[1][2] == 0: # if the direction is has a diagonal direction after it, pop them both out, as the path after the diagonal direction now won't be complete
            arrows_matrix[paths[0][0], paths[0][1], paths[0][2]] = False
            arrows_matrix[paths[1][0], paths[1][1], paths[1][2]] = False
            paths.pop(0)
            paths.pop(0)

    return possible_layouts


#  DNA sequence alignment using Needleman-Wunsch Algorithm
def global_alignment(dir1, dir2, match_score, mismatch_penalty, gap_penalty):
    # Read sequences from files
    seq1, seq2 = parse_input(dir1, dir2)

    # Fill score and arrows matrices with initial values
    initialization_matrix, arrows_matrix = initialize(seq1, seq2, gap_penalty)

    # Calculate scores and directions of the remaining positions in the matrices
    score_matrix, arrows_matrix = matrix_fill(initialization_matrix, arrows_matrix,
                                              match_score, mismatch_penalty, gap_penalty, seq1, seq2)

    # Find all possible layouts that results in the best score
    possible_layouts = find_paths(score_matrix, arrows_matrix, seq1, seq2)

    print("Best score is {}\n".format(score_matrix[len(seq2), len(seq1)]))
    for layout in possible_layouts:
        print("Possible layout: \n{}\n{}\n".format(layout[0], layout[1]))


# Test case for the 2 sequences in lecture 7 with the scores from assignment 1
global_alignment('sequence1.txt', 'sequence2.txt', 4, -1, -3)

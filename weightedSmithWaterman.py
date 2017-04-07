#!/bin/python
# (c) 2013 Ryan Boehning
# (c) 2017 Daniel McSkimming


'''A Python implementation of the Smith-Waterman algorithm for local alignment
of nucleotide sequences weighted by base quality (of seq2).
'''


def do_one(args):
    try:
        parse_cmd_line()
    except ValueError as err:
        print('error:', err)
        return

    # The scoring matrix contains an extra row and column for the gap (-), hence
    # the +1 here.
    rows = len(args.seq1) + 1
    cols = len(args.seq2) + 1

    # Initialize the scoring matrix.
    score_matrix, start_pos = create_score_matrix(rows, cols, args.match,\
            args.mismatch, args.gap)

    # Traceback. Find the optimal path through the scoring matrix. This path
    # corresponds to the optimal local sequence alignment.
    seq1_aligned, seq2_aligned = traceback(score_matrix, start_pos, args.seq1, args.seq2)
    assert len(seq1_aligned) == len(seq2_aligned), 'aligned strings are not the same size'

    # Pretty print the results. The printing follows the format of BLAST results
    # as closely as possible.
    alignment_str, idents, gaps, mismatches = alignment_string(seq1_aligned, seq2_aligned)
    alength = len(seq1_aligned)
    print()
    print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents,
          alength, idents / alength, gaps, alength, gaps / alength))
    print()
    for i in range(0, alength, 60):
        seq1_slice = seq1_aligned[i:i+60]
        print('Query  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
        print('             {0}'.format(alignment_str[i:i+60]))
        seq2_slice = seq2_aligned[i:i+60]
        print('Sbjct  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
        print()


def create_score_matrix(seq1, seq2, match=2, mismatch=-1, gap=-1):
    # pass in seq1 and seq2, not rows,cols
    '''Create a matrix of scores representing trial alignments of the two sequences.

    Sequence alignment can be treated as a graph search problem. This function
    creates a graph (2D matrix) of scores, which are based on trial alignments
    of different base pairs. The path with the highest cummulative score is the
    best alignment.
    '''
    clen, rlen = len(seq1)+1, len(seq2)+1
    score_matrix = [[0 for col in range(clen)] for row in range(rlen)]
    # Fill the scoring matrix.
    max_score, score = 0, 0
    seq2_qual = 1.0
    max_pos   = None    # The row and columbn of the highest score in matrix.
    for i in range(1, rlen):
        for j in range(1, clen):
            xprev_y, x_yprev, xprev_yprev = score_matrix[i-1][j], score_matrix[i][j-1], score_matrix[i-1][j-1]
            score = calc_score(xprev_y, x_yprev, xprev_yprev, seq1[i-1], seq2[j-1], seq2_qual, match, mismatch, gap) 
            if score > max_score:
                max_score = score
                max_pos   = (i, j)
            score_matrix[i][j] = score

    assert max_pos is not None, 'the x, y position with the highest score was not found'

    return score_matrix, max_pos


def calc_score(xpy, xyp, xpyp, a, b, b_qual, match, mismatch, gap):
    '''Calculate score for a given x, y position in the scoring matrix.

    The score is based on the up, left, and upper-left neighbors.
    '''
    #need to scale match/mismatch by base quality
    similarity = b_qual*match if a == b else b_qual*mismatch

    diag_score = xpyp + similarity
    up_score   = xpy + gap
    left_score = xyp + gap

    return max(0, diag_score, up_score, left_score)


def traceback(score_matrix, start_pos, seq1, seq2):
    '''Find the optimal path through the matrix.

    This function traces a path from the bottom-right to the top-left corner of
    the scoring matrix. Each move corresponds to a match, mismatch, or gap in one
    or both of the sequences being aligned. Moves are determined by the score of
    three adjacent squares: the upper square, the left square, and the diagonal
    upper-left square.

    WHAT EACH MOVE REPRESENTS
        diagonal: match/mismatch
        up:       gap in sequence 1
        left:     gap in sequence 2
    '''

    END, DIAG, UP, LEFT = range(4)
    aligned_seq1 = []
    aligned_seq2 = []
    x, y         = start_pos
    xpyp, xpy, xyp = score_matrix[x-1][y-1], score_matrix[x-1][y],\
            score_matrix[x][y-1]
    move = next_move(xpyp, xpy, xyp)
    while move != END:
        if move == DIAG:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append(seq2[y - 1])
            x -= 1
            y -= 1
        elif move == UP:
            aligned_seq1.append(seq1[x - 1])
            aligned_seq2.append('-')
            x -= 1
        else:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[y - 1])
            y -= 1
        xpyp, xpy, xyp = score_matrix[x-1][y-1], score_matrix[x-1][y],\
                score_matrix[x][y-1]
        move = next_move(xpyp, xpy, xyp)

    aligned_seq1.append(seq1[x - 1])
    aligned_seq2.append(seq2[y - 1])

    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


def next_move(diag, up, left):
    if diag >= up and diag >= left:     # Tie goes to the DIAG move.
        return 1 if diag != 0 else 0    # 1 signals a DIAG move. 0 signals the end.
    elif up > diag and up >= left:      # Tie goes to UP move.
        return 2 if up != 0 else 0      # UP move or end.
    elif left > diag and left > up:
        return 3 if left != 0 else 0    # LEFT move or end.
    else:
        # Execution should not reach here.
        raise ValueError('invalid move during traceback')
    return


def alignment_string(aligned_seq1, aligned_seq2):
    '''Construct a special string showing identities, gaps, and mismatches.

    This string is printed between the two aligned sequences and shows the
    identities (|), gaps (-), and mismatches (:). As the string is constructed,
    it also counts number of identities, gaps, and mismatches and returns the
    counts along with the alignment string.

    AAGGATGCCTCAAATCGATCT-TTTTCTTGG-
    ::||::::::||:|::::::: |:  :||:|   <-- alignment string
    CTGGTACTTGCAGAGAAGGGGGTA--ATTTGG
    '''
    # Build the string as a list of characters to avoid costly string
    # concatenation.
    idents, gaps, mismatches = 0, 0, 0
    alignment_string = []
    for base1, base2 in zip(aligned_seq1, aligned_seq2):
        if base1 == base2:
            alignment_string.append('|')
            idents += 1
        elif '-' in (base1, base2):
            alignment_string.append(' ')
            gaps += 1
        else:
            alignment_string.append(':')
            mismatches += 1

    return ''.join(alignment_string), idents, gaps, mismatches


def print_matrix(matrix):
    '''Print the scoring matrix.

    ex:
    0   0   0   0   0   0
    0   2   1   2   1   2
    0   1   1   1   1   1
    0   0   3   2   3   2
    0   2   2   5   4   5
    0   1   4   4   7   6
    '''
    for row in matrix:
        print('\t'.join('{0:>6}'.format(x)for x in row))
        #for col in row:
            #print('{0:>6}'.format(col))
    return

if __name__ == '__main__':
    import argparse
    #use argparse to take sequences, validate format, score and display alignment

import unittest
import weightedSmithWaterman as wSW

class ScoreMatrixTest(unittest.TestCase):
    '''Compare the matrix produced by create_score_matrix() with a known matrix.'''
    def test_nonweighted_matrix(self):
        # From Wikipedia (en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
        #                -   A   C   A   C   A   C   T   A
        known_matrix = [[0,  0,  0,  0,  0,  0,  0,  0,  0],  # -
                        [0,  2,  1,  2,  1,  2,  1,  0,  2],  # A
                        [0,  1,  1,  1,  1,  1,  1,  0,  1],  # G
                        [0,  0,  3,  2,  3,  2,  3,  2,  1],  # C
                        [0,  2,  2,  5,  4,  5,  4,  3,  4],  # A
                        [0,  1,  4,  4,  7,  6,  7,  6,  5],  # C
                        [0,  2,  3,  6,  6,  9,  8,  7,  8],  # A
                        [0,  1,  4,  5,  8,  8, 11, 10,  9],  # C
                        [0,  2,  3,  6,  7, 10, 10, 10, 12]]  # A

        seq1 = 'AGCACACA'
        seq2 = 'ACACACTA'

        matrix_to_test, max_pos = wSW.create_score_matrix(seq1, seq2, 2, -1, -1)
        wSW.print_matrix(matrix_to_test)
        seq1_aligned, seq2_aligned = wSW.traceback(matrix_to_test, max_pos, seq1, seq2)
        self.assertEqual(len(seq1_aligned), len(seq2_aligned))
        alignment_str, idents, gaps, mismatches = wSW.alignment_string(seq1_aligned, seq2_aligned)
        alength = len(seq1_aligned)
        print('\n')
        print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents,
            alength, float(idents) / alength, gaps, alength, float(gaps) / alength))
        print('\n')
        for i in range(0, alength, 60):
            seq1_slice = seq1_aligned[i:i+60]
            print('Query  {0:<4}  {1}  {2:<4}'.format(i + 1, seq1_slice, i + len(seq1_slice)))
            print('             {0}'.format(alignment_str[i:i+60]))
            seq2_slice = seq2_aligned[i:i+60]
            print('Sbjct  {0:<4}  {1}  {2:<4}'.format(i + 1, seq2_slice, i + len(seq2_slice)))
            print('\n')
            self.assertEqual(known_matrix, matrix_to_test)
        return

if __name__ == '__main__':
    unittest.main()

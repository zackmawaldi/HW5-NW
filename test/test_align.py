# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    Unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """

    # MYQR
    seq1, _ = read_fasta("./data/test_seq1.fa")

    # MQR
    seq2, _ = read_fasta("./data/test_seq2.fa")

    true_align = np.array([[  0, -11, -12, -13],
                           [-11,   5,  -6,  -7],
                           [-12,  -6,   4,  -7],
                           [-13,  -7,  -1,   5],
                           [-14,  -8,  -6,   4]])
    
    true_is_gap = np.array([[0, 1, 1, 1],
                            [1, 0, 1, 1],
                            [1, 1, 0, 1],
                            [1, 1, 0, 0],
                            [1, 1, 0, 0]])
    
    true_where_from = np.array([[None, None,  None,   None],
                            [None, (0, 0), (1, 1), (1, 2)],
                            [None, (1, 1), (1, 1), (2, 2)],
                            [None, (2, 1), (2, 1), (2, 2)],
                            [None, (3, 1), (3, 1), (3, 2)]], dtype=object)
    

    nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10, -1)
    nw.align(seq1, seq2)

    assert np.array_equal( nw._align_matrix , true_align )
    assert np.array_equal( nw._is_gap , true_is_gap )
    assert np.array_equal( nw._where_from , true_where_from )
    

def test_nw_backtrace():
    """
    Unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """

    # MAVHQLIRRP
    seq3, _ = read_fasta("./data/test_seq3.fa")

    # MQLIRHP
    seq4, _ = read_fasta("./data/test_seq4.fa")

    nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', -10, -1)
    score, seqA_align, seqB_align = nw.align(seq3, seq4)

    assert round(score) == 17
    assert seqA_align   == 'MAVHQLIRRP'
    assert seqB_align   == 'M---QLIRHP'
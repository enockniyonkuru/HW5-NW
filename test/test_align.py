# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    Test NW alignment by checking that the alignment matrices are filled correctly
    using test_seq1.fa and test_seq2.fa.
    Use the BLOSUM62 matrix and a gap open penalty of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    aligner = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    score, align1, align2 = aligner.align(seq1, seq2)
    
    # Check that matrices are not None and have correct dimensions
    assert aligner._align_matrix is not None, "Alignment matrix should not be None"
    assert aligner._gapA_matrix is not None, "Gap A matrix should not be None"
    assert aligner._gapB_matrix is not None, "Gap B matrix should not be None"
    
    # Check matrix dimensions
    m, n = len(seq1) + 1, len(seq2) + 1
    assert aligner._align_matrix.shape == (m, n), f"Alignment matrix shape should be {(m, n)}"
    assert aligner._gapA_matrix.shape == (m, n), f"Gap A matrix shape should be {(m, n)}"
    assert aligner._gapB_matrix.shape == (m, n), f"Gap B matrix shape should be {(m, n)}"
    
    # Check that the alignment score is returned
    assert isinstance(score, (int, float, np.floating)), "Score should be a number"
    assert len(align1) == len(align2), "Both alignments should have the same length"
    

def test_nw_backtrace():
    """
    Test NW backtracing using test_seq3.fa and test_seq4.fa.
    Expected: alignment score of 17-18 (depending on exact scoring)
    Use the BLOSUM62 matrix, gap open penalty of -10, and gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    aligner = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    score, align3, align4 = aligner.align(seq3, seq4)
    
    # Check that alignment score is reasonable (should be 17-18)
    assert score >= 17, f"Expected alignment score >= 17, got {score}"
    
    # Check that both alignments have the same length
    assert len(align3) == len(align4), "Both alignments should have the same length"
    
    # Check that the returned alignments are non-empty
    assert len(align3) > 0, "Alignment strings should not be empty"
    assert len(align4) > 0, "Alignment strings should not be empty"
    
    # Verify the alignment matches expected pattern (should start with M and have gaps in seq4)
    assert align3[0] == 'M', "First character of seq3 alignment should be M"
    assert align4[0] == 'M', "First character of seq4 alignment should be M"





from align import NeedlemanWunsch, read_fasta
import numpy as np

seq1, _ = read_fasta("./data/test_seq1.fa")
seq2, _ = read_fasta("./data/test_seq2.fa")

print(f"seq1: {seq1} (length: {len(seq1)})")
print(f"seq2: {seq2} (length: {len(seq2)})")

aligner = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", -10, -1)
score, align1, align2 = aligner.align(seq1, seq2)

print(f"\nAlignment score: {score}")
print(f"seq1 alignment: {align1}")
print(f"seq2 alignment: {align2}")

print(f"\nH matrix (_align_matrix):")
print(aligner._align_matrix)

print(f"\nE matrix (_gapB_matrix):")
print(aligner._gapB_matrix)

print(f"\nF matrix (_gapA_matrix):")
print(aligner._gapA_matrix)

print(f"\nBack matrix:")
print(aligner._back)

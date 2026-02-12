from align import NeedlemanWunsch, read_fasta
import numpy as np

seq3, _ = read_fasta("./data/test_seq3.fa")
seq4, _ = read_fasta("./data/test_seq4.fa")

aligner = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)

# Check BLOSUM62 scores
print("BLOSUM62 scores:")
print(f"M-M: {aligner.sub_dict.get(('M', 'M'))}")
print(f"A-A: {aligner.sub_dict.get(('A', 'A'))}")
print(f"Q-Q: {aligner.sub_dict.get(('Q', 'Q'))}")
print(f"L-L: {aligner.sub_dict.get(('L', 'L'))}")
print(f"I-I: {aligner.sub_dict.get(('I', 'I'))}")
print(f"R-R: {aligner.sub_dict.get(('R', 'R'))}")
print(f"R-H: {aligner.sub_dict.get(('R', 'H'))}")
print(f"P-P: {aligner.sub_dict.get(('P', 'P'))}")

score, align3, align4 = aligner.align(seq3, seq4)

print(f"\nAlignment score: {score}")
print(f"seq3: {align3}")
print(f"seq4: {align4}")

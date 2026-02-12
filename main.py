# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    # Read sequences
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/Tursiops_truncatus_BRD2.fa")
    
    # Create aligner with BLOSUM62 matrix, gap_open=-10, gap_extend=-1
    aligner = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    
    # List of species with their sequences and headers
    species_data = [
        ("Gallus gallus", gg_seq),
        ("Mus musculus", mm_seq),
        ("Balaeniceps rex", br_seq),
        ("Tursiops truncatus", tt_seq)
    ]
    
    # Align each species to human and store results
    alignments = []
    for species_name, species_seq in species_data:
        score, seqA_align, seqB_align = aligner.align(hs_seq, species_seq)
        alignments.append((species_name, score))
    
    # Sort by score (descending - most similar first)
    alignments.sort(key=lambda x: x[1], reverse=True)
    
    # Print species ordered by similarity to human
    print("Species ordered by similarity to Homo sapiens BRD2 (most similar to least similar):")
    for i, (species_name, score) in enumerate(alignments, 1):
        print(f"{i}. {species_name}")
    
    print("\nAlignment scores for each species aligned to human BRD2:")
    for species_name, score in alignments:
        print(f"{species_name}: {score}")
    

if __name__ == "__main__":
    main()

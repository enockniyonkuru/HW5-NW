# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm with affine gap penalties.
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # Get sequence lengths
        m = len(seqA)
        n = len(seqB)
        
        # Initialize matrix private attributes for use in alignment
        # H matrix: main alignment scores
        self._align_matrix = np.full((m + 1, n + 1), -np.inf)
        # E matrix: gap extension for seqB (horizontal)
        self._gapB_matrix = np.full((m + 1, n + 1), -np.inf)
        # F matrix: gap extension for seqA (vertical)
        self._gapA_matrix = np.full((m + 1, n + 1), -np.inf)
        
        # Backtracing matrices to track which matrix produced the max score
        self._back = np.zeros((m + 1, n + 1), dtype=int)  # 0=H, 1=E, 2=F
        self._back_A = np.zeros((m + 1, n + 1), dtype=int)  # for gapA
        self._back_B = np.zeros((m + 1, n + 1), dtype=int)  # for gapB
        
        # Initialize first row and column
        self._align_matrix[0, 0] = 0
        for i in range(1, m + 1):
            self._gapA_matrix[i, 0] = self.gap_open + (i - 1) * self.gap_extend
            self._align_matrix[i, 0] = self._gapA_matrix[i, 0]
            self._back[i, 0] = 2  # Coming from F matrix
            
        for j in range(1, n + 1):
            self._gapB_matrix[0, j] = self.gap_open + (j - 1) * self.gap_extend
            self._align_matrix[0, j] = self._gapB_matrix[0, j]
            self._back[0, j] = 1  # Coming from E matrix
        
        # Fill in the dynamic programming matrices
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                # Get the substitution score for the current residues
                sub_score = self.sub_dict.get((seqA[i-1], seqB[j-1]), -np.inf)
                
                # Calculate E matrix (gap in seqA, gap extends horizontally in seqB)
                e1 = self._align_matrix[i, j-1] + self.gap_open
                e2 = self._gapB_matrix[i, j-1] + self.gap_extend
                if e1 >= e2:
                    self._gapB_matrix[i, j] = e1
                    self._back_B[i, j] = 0  # From H
                else:
                    self._gapB_matrix[i, j] = e2
                    self._back_B[i, j] = 1  # From E
                
                # Calculate F matrix (gap in seqB, gap extends vertically in seqA)
                f1 = self._align_matrix[i-1, j] + self.gap_open
                f2 = self._gapA_matrix[i-1, j] + self.gap_extend
                if f1 >= f2:
                    self._gapA_matrix[i, j] = f1
                    self._back_A[i, j] = 0  # From H
                else:
                    self._gapA_matrix[i, j] = f2
                    self._back_A[i, j] = 2  # From F
                
                # Calculate H matrix (alignment score)
                h1 = self._align_matrix[i-1, j-1] + sub_score
                h2 = self._gapB_matrix[i, j]
                h3 = self._gapA_matrix[i, j]
                
                # Find maximum value
                max_val = max(h1, h2, h3)
                self._align_matrix[i, j] = max_val
                
                # Record which matrix produced the max
                if max_val == h1:
                    self._back[i, j] = 0  # From H (match/mismatch)
                elif max_val == h2:
                    self._back[i, j] = 1  # From E (gap in A)
                else:
                    self._back[i, j] = 2  # From F (gap in B)
        
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Store the alignment score
        self.alignment_score = self._align_matrix[len(self._seqA), len(self._seqB)]
        
        # Initialize alignment strings
        seqA_align = ""
        seqB_align = ""
        
        # Start from bottom-right corner
        i = len(self._seqA)
        j = len(self._seqB)
        
        # Current matrix: 0=H, 1=E, 2=F
        current_matrix = 0
        
        # Backtrace through the matrices
        while i > 0 or j > 0:
            if current_matrix == 0:  # In H matrix
                if i == 0:
                    # Top edge: must be gap in seqA
                    seqA_align = "-" + seqA_align
                    seqB_align = self._seqB[j-1] + seqB_align
                    j -= 1
                    current_matrix = 1  # Move to E matrix
                elif j == 0:
                    # Left edge: must be gap in seqB
                    seqA_align = self._seqA[i-1] + seqA_align
                    seqB_align = "-" + seqB_align
                    i -= 1
                    current_matrix = 2  # Move to F matrix
                else:
                    # Determine which matrix gave the maximum
                    back_val = self._back[i, j]
                    
                    if back_val == 0:  # Came from diagonal (match/mismatch)
                        seqA_align = self._seqA[i-1] + seqA_align
                        seqB_align = self._seqB[j-1] + seqB_align
                        i -= 1
                        j -= 1
                    elif back_val == 1:  # Came from E matrix (gap in A)
                        current_matrix = 1
                    else:  # Came from F matrix (gap in B)
                        current_matrix = 2
                        
            elif current_matrix == 1:  # In E matrix (gap in A)
                seqA_align = "-" + seqA_align
                seqB_align = self._seqB[j-1] + seqB_align
                j -= 1
                
                # Check where to go next in E matrix
                if j == 0:
                    current_matrix = 0
                else:
                    back_val = self._back_B[i, j]
                    if back_val == 0:  # Go back to H
                        current_matrix = 0
                    else:  # Stay in E
                        current_matrix = 1
                        
            else:  # In F matrix (gap in B)
                seqA_align = self._seqA[i-1] + seqA_align
                seqB_align = "-" + seqB_align
                i -= 1
                
                # Check where to go next in F matrix
                if i == 0:
                    current_matrix = 0
                else:
                    back_val = self._back_A[i, j]
                    if back_val == 0:  # Go back to H
                        current_matrix = 0
                    else:  # Stay in F
                        current_matrix = 2
        
        self.seqA_align = seqA_align
        self.seqB_align = seqB_align
        
        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header

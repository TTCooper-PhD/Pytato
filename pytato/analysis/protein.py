import numpy as np
import re

class ProteinAnalysis:
    """
    Methods for protein-level analysis
    """
    @staticmethod
    def Pep2Pro(protein, peptides):
        """
        Calculate protein coverage from peptide matches
        
        Parameters:
        protein (str): Protein sequence
        peptides (list): List of peptide sequences
        
        Returns:
        float: Coverage ratio (0-1)
        """
        protein = re.sub(r'[^A-Z]', '', protein)
        mask = np.zeros(len(protein), dtype=np.int8)
        for peptide in peptides:
            indices = [m.start() for m in re.finditer(
                '(?={})'.format(re.sub(r'[^A-Z]', '', peptide)), protein)]
            for i in indices:
                mask[i:i + len(peptide)] = 1
        return mask.sum() / mask.size  # Fixed typo here 
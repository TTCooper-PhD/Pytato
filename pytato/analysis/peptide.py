import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import ttest_ind
import re
from collections import deque
from Bio import SeqIO
import requests
from io import StringIO
import gzip
import itertools as it
from pytato.core.scales import Scales



class Peptide:
    def __init__(self):
        self.pH=2.0
        self.missed = 0
        self.min_len = 6
        self.max_len = 100
        self.enzyme = "trypsin"
        self.rules = { 'arg-c':         r'R',
        'asp-n':         r'\w(?=D)',
        'bnps-skatole' : r'W',
        'caspase 1':     r'(?<=[FWYL]\w[HAT])D(?=[^PEDQKR])',
        'caspase 2':     r'(?<=DVA)D(?=[^PEDQKR])',
        'caspase 3':     r'(?<=DMQ)D(?=[^PEDQKR])',
        'caspase 4':     r'(?<=LEV)D(?=[^PEDQKR])',
        'caspase 5':     r'(?<=[LW]EH)D',
        'caspase 6':     r'(?<=VE[HI])D(?=[^PEDQKR])',
        'caspase 7':     r'(?<=DEV)D(?=[^PEDQKR])',
        'caspase 8':     r'(?<=[IL]ET)D(?=[^PEDQKR])',
        'caspase 9':     r'(?<=LEH)D',
        'caspase 10':    r'(?<=IEA)D',
        'chymotrypsin high specificity' : r'([FY](?=[^P]))|(W(?=[^MP]))',
        'chymotrypsin low specificity':
            r'([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))',
        'clostripain':   r'R',
        'cnbr':          r'M',
        'enterokinase':  r'(?<=[DE]{3})K',
        'factor xa':     r'(?<=[AFGILTVM][DE]G)R',
        'formic acid':   r'D',
        'glutamyl endopeptidase': r'E',
        'granzyme b':    r'(?<=IEP)D',
        'hydroxylamine': r'N(?=G)',
        'iodosobenzoic acid': r'W',
        'lysc':          r'K',
        'ntcb':          r'\w(?=C)',
        'pepsin ph1.3':  r'((?<=[^HKR][^P])[^R](?=[FL][^P]))|'
                        r'((?<=[^HKR][^P])[FL](?=\w[^P]))',
        'pepsin ph2.0':  r'((?<=[^HKR][^P])[^R](?=[FLWY][^P]))|'
                        r'((?<=[^HKR][^P])[FLWY](?=\w[^P]))',
        'proline endopeptidase': r'(?<=[HKR])P(?=[^P])',
        'proteinase k':  r'[AEFILTVWY]',
        'staphylococcal peptidase i': r'(?<=[^E])E',
        'thermolysin':   r'[^DE](?=[AFILMV])',
        'thrombin':      r'((?<=G)R(?=G))|'
                        r'((?<=[AFGILTVM][AFGILTVWA]P)R(?=[^DE][^DE]))',
        'trypsin':       r'([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))',
        'trypsin_exception': r'((?<=[CD])K(?=D))|((?<=C)K(?=[HY]))|((?<=C)R(?=K))|((?<=R)R(?=[HR]))'}

    @staticmethod
    def modify_peptides(peptides_df, amino_acid='C', mass_shift=105.0578, charge_shift=1):
        """Modify peptides containing specific amino acid by adding mass and charge."""
        modified_df = peptides_df.copy()
        

        # Count occurrences of specified amino acid
        modified_df[f'{amino_acid}_count'] = modified_df['peptide'].str.count(amino_acid)
        
        # Add mass shift and charge for each occurrence
        modified_df['m/z'] = modified_df.apply(
            lambda x: x['m/z'] + (mass_shift * x[f'{amino_acid}_count']) 
            if pd.notnull(x['m/z']) else x['m/z'], 
            axis=1
        )
        modified_df['z'] = modified_df['z'] + (modified_df[f'{amino_acid}_count'] * charge_shift)
        
        # Drop the temporary count column
        modified_df.drop(f'{amino_acid}_count', axis=1, inplace=True)
        
        return modified_df
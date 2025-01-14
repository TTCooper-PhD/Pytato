from Bio import SeqIO
import requests
import gzip
from io import StringIO
import pandas as pd
import re

class Silico:
    def __init__(self):
        self.pH = 2.0
        self.missed = 0
        self.min_len = 6
        self.max_len = 100
        self.enzyme = "trypsin"
        self.rules = {...}  # Move enzyme rules here
        
    def import_fasta(self, source=None, local_path=None):
        # Move existing import_fasta method here
        pass

    def cleave_sequence(self, sequence, enzyme=None, exception=None):
        # Move existing cleave_sequence method here
        pass 
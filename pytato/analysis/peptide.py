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



class Silico:
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

    def update_rules(self, new_enzyme=None, new_rule=None):
        """
        Update the enzyme rules dictionary with a new enzyme and its corresponding cleavage rule. Rules must be regex.
        """
        if new_enzyme and new_rule:
            self.rules[new_enzyme] = new_rule

    def import_fasta(self, source=None, local_path=None, default_url="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"):
        """
        Import a FASTA file from UniProt or a local file and convert it into a DataFrame.

        Parameters:
        - source (str): URL of the FASTA file. If None, defaults to a human proteome file from UniProt.
        - local_path (str): Path to a local FASTA file. If specified, `source` is ignored.

        Returns:
        - pd.DataFrame: DataFrame with columns ['UniprotID', 'Gene', 'Peptide']
        """
        ingredients = []

        if local_path:
            print("Importing FASTA file from local path...")
            with open(local_path, 'rt') as file:
                for record in SeqIO.parse(file, "fasta"):
                    ingredients.append([record.id, str(record.seq)])
        else:
            print("Downloading and importing FASTA file from URL...")
            source = source if source else default_url
            response = requests.get(source)
            fasta_text = gzip.decompress(response.content).decode('utf-8')
            for record in SeqIO.parse(StringIO(fasta_text), "fasta"):
                ingredients.append([record.id, str(record.seq)])

        # Convert to DataFrame and parse details
        recipie = pd.DataFrame(ingredients, columns=['ID', 'Peptide'])
        recipie[['db', 'UniprotID', 'ID2']] = recipie['ID'].str.split('|', expand=True)
        recipie[['Gene', 'Identification']] = recipie['ID2'].str.split('_', expand=True)
        recipie.drop(columns=['ID', 'ID2', 'db'], inplace=True)

        print("FASTA import complete!")
        return recipie[['UniprotID', 'Gene', 'Peptide']]

    def cleave_sequence(self, sequence, enzyme=None, exception=None):
        """
        Simulate enzymatic digestion of a protein sequence based on specified rules and parameters.
        
        Parameters:
        - sequence (str): The protein sequence to be digested.
        - enzyme (str): The enzyme to use for digestion. If None, uses the class default enzyme.
        - exception (str): A rule for exceptions in the cleavage pattern, if any.

        Returns:
        - list[str]: A list of resulting peptide sequences.
        """
        peptides = []
        rule = self.rules.get(enzyme, self.rules.get(self.enzyme))  # Use specified enzyme or class default
        if rule is None:
            raise ValueError("Enzyme & Rule Unknown. Please use update_rules to add new enzyme rules.")
        
        exception_rule = self.rules.get(exception, None) if exception else None
        ml = self.missed + 2
        trange = range(ml)  # returns range of 0 to ml-1
        cleavage_sites = deque([0], maxlen=ml)  # deque to store cleavage sites
        
        if exception_rule:
            exceptions = {x.end() for x in re.finditer(exception_rule, sequence)}
        
        for i in it.chain([x.end() for x in re.finditer(rule, sequence)], [None]):
            if exception_rule and i in exceptions:
                continue  # Skip exceptions
            cleavage_sites.append(i)
            for j in trange[:len(cleavage_sites) - 1]:
                seq = sequence[cleavage_sites[j]:cleavage_sites[-1]]
                if self.min_len <= len(seq) <= self.max_len:
                    peptides.append(seq)
        
        return peptides

    def clean_data(self, df, filters=[], t_id=None, t_value=0, acid=["J", "Z"], labels=[]):
        """
        Clean and preprocess peptide data from a DataFrame with flexible filtering.

        Parameters:
        - df (pd.DataFrame): The DataFrame containing peptide data.
        - filters (list): A list of strings representing column name patterns to filter the DataFrame in sequential order.
        - t_id (str): A regex pattern to identify columns to check against t_value.
        - t_value (float): The threshold value for filtering columns identified by t_id.
        - acid (list): Amino acids to exclude.
        - labels (list): Column names to preserve during filtering.

        Returns:
        - pd.DataFrame: The processed DataFrame after applying all filters and conditions.
        """
        tag = df[labels]

        # Sequentially apply filters
        filtered_df = df.copy()
        for f in filters:
            if f != "Pass":
                filtered_df = filtered_df.filter(like=str(f), axis=1)
                filtered_df = pd.concat([filtered_df, tag], axis=1)

        # Clean peptide sequences
        filtered_df["Peptide"] = filtered_df["Peptide"].str.replace('\W+|\d+', "")
        filtered_df["Peptide"] = filtered_df["Peptide"].apply(lambda x: x.strip("[]"))

        if t_id:
            # Filter columns based on t_id and t_value
            blade = filtered_df.filter(regex=t_id, axis='columns')
            trim = blade <= t_value
            trimmings = trim.all(axis=1)
            filtered_df = filtered_df.loc[~trimmings]

        filtered_df.reset_index(drop=True, inplace=True)

        return filtered_df

    def generate_samples(self, df, target="Peptide", identifier="Gene", enzyme=None, min_length=7, exception=None, max_length=100, pH=None, min_charge=2.0):
        """
        Generate artificial datasets based on enzymatic digestion rules and analyze peptide properties.

        Parameters:
        - df (pd.DataFrame): DataFrame containing the target peptides and identifiers.
        - target (str): Column name for peptide sequences.
        - identifier (str): Column name for the peptide identifiers (e.g., gene names).
        - enzyme (str): Enzyme used for digestion. Uses the class's default enzyme if None.
        - min_length (int): Minimum length of peptides to consider.
        - exception (str): Exception rule for the enzyme cleavage.
        - max_length (int): Maximum length of peptides to consider.
        - pH (float): pH value to calculate peptide charge.
        - min_charge (float): Minimum charge threshold for peptides to be included in the results.

        Returns:
        - list[dict]: A list of dictionaries, each representing a peptide with its properties.
        """
        peptide_properties = []
        enzyme = enzyme if enzyme else self.enzyme
        pH = pH if pH else self.pH
        print(f'Processing your order with {enzyme}-cut proteins.')

        for index, row in df.iterrows():
            gene = row[identifier]
            sequence = row[target]
            peptides = self.cleave_sequence(sequence, enzyme=enzyme, exception=exception)

            for peptide in peptides:
                if min_length <= len(peptide) <= max_length:
                    properties = {
                        'gene': gene,
                        'peptide': peptide,
                        'Length': len(peptide),
                        'aa_comp': Scales.peptide_inspector(peptide), 
                        'neutral_z': Scales.z_neutral_ph(peptide),
                        'z':Scales.calculate_peptide_charge(peptide,pH), 
                        'Mass': Scales.calculate_mass(peptide,), 
                        'GRAVY':Scales.peptide_gravy(peptide)
                    }

                    if properties["z"] >= min_charge:
                        if properties["z"] > 0:
                            properties.update({'m/z': properties["Mass"] / properties['z']})
                        peptide_properties.append(properties)

        print(f'Generated {len(peptide_properties)} peptides matching criteria.')
        return peptide_properties

    def capture_flanking_sequences(self, protein_sequence, peptide_sequence, flank_length):
        """
        Capture the amino acid sequences flanking the cleavage site of a given peptide within a protein.

        Parameters:
        protein_sequence (str): The full amino acid sequence of the protein.
        peptide_sequence (str): The amino acid sequence of the identified peptide.
        flank_length (int): The number of amino acids to capture on each side of the cleavage site.

        Returns:
        dict: A dictionary with two keys, 'n_cut' and 'c_cut', containing sequences flanking the cleavage site.
        """
        peptide_start = protein_sequence.find(peptide_sequence)
        if peptide_start == -1:
            raise ValueError("Peptide sequence not found in protein sequence")

        # Find all cleavage points in the peptide sequence
        cleavage_points = [m.start() for m in re.finditer(self.rule, peptide_sequence)]
        if self.exception:
            cleavage_points = [cp for cp in cleavage_points if not re.match(self.exception, peptide_sequence[cp:])]

        flanking_sequences = []

        for cp in cleavage_points:
            # Adjust the cleavage point to the full protein sequence
            adjusted_cp = peptide_start + cp

            # Capture flanking sequences
            n_cut_start = max(0, adjusted_cp - flank_length)
            c_cut_end = min(len(protein_sequence), adjusted_cp + flank_length + 1)

            n_cut = protein_sequence[n_cut_start:adjusted_cp]
            c_cut = protein_sequence[adjusted_cp + 1:c_cut_end]

            flanking_sequences.append({'n_cut': n_cut, 'c_cut': c_cut})

        return flanking_sequences
    
    def Pep2Pro(self, protein, peptides):
        protein = re.sub(r'[^A-Z]', '', protein)
        mask = np.zeros(len(protein), dtype=np.int8)
        for peptide in peptides:
            indices = [m.start() for m in re.finditer(
                '(?={})'.format(re.sub(r'[^A-Z]', '', peptide)), protein)]
            for i in indices:
                mask[i:i + len(peptide)] = 1
        return mask.sum / mask.size


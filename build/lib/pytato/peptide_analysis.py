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




class Silico:
    def __init__(self):
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

    def generate_samples(self, df, target="Peptide", identifier="Gene", enzyme=None, min_length=7, exception=None, max_length=100, pH=2.0, min_charge=2.0):
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
                        'z': Scales.z_neutral_ph(peptide),  
                        'Mass': Scales.calculate_mass(peptide), 
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
    

class Scales:
    # Base masses for the amino acids
    base_mass = {  
        "A": 71.037114,
        "R": 156.101111,
        "N": 114.042927,
        "D": 115.026943,
        "C": 103.009185,
        "Q": 129.042593,
        "E": 128.058578,
        "G": 57.021464,
        "H": 137.058912,
        "I": 113.084064,
        "L": 113.084064,
        "K": 128.094963,
        "M": 131.040485,
        "F": 147.068414,
        "P": 97.052764,
        "S": 87.032028,
        "T": 101.047679,
        "W": 186.079313,
        "Y": 163.06332,
        "V": 99.068414,
    }

    @staticmethod
    def parse_modifications(peptide: str) -> dict:
        """
        Parse modifications in a peptide sequence and return a dictionary with positions and modified masses.

        Parameters:
        peptide (str): The peptide sequence with modifications (e.g., "C[+57.02146]").

        Returns:
        dict: A dictionary with positions as keys and their respective masses (base + modification) as values.

        # Example usage:
        scales = Scales()
        modifications = scales.parse_modifications("LVC[+57.02146]TALQW")
        print(modifications)

        Returns:
        {0: 113.084064, 1: 99.068414, '2_C': 160.030645, 3: 101.047679, 4: 71.037114, 5: 113.084064, 6: 129.042593, 7: 186.079313}

        ## Note that 2_C is a non-numeric (unique) key, giving the position of the amino acid (2) and the amino acid itself (C)
        """
        # Regular expression to match modifications (e.g., "C[+57.02146]")
        pattern = re.compile(r'([A-Z])(\[\+?(-?\d+\.\d+)\])?')
        mods = {}
        pos = 0

        for match in pattern.finditer(peptide):
            aa, mod_str, mod_mass = match.groups()
            if mod_mass:
                # Use the position in the peptide and the amino acid to create a unique key
                key = f"{pos}_{aa}"
                mods[key] = Scales.base_mass.get(aa, 0.0) + float(mod_mass)
            else:
                # If no modification, just add the base mass
                mods[pos] = Scales.base_mass.get(aa, 0.0)
            pos += 1
        
        return mods

    @staticmethod
    def calculate_mass(peptide: str) -> float:
        """
        Calculate the total monoisotopic mass of a peptide sequence, including any modifications.

        Parameters:
        peptide (str): The amino acid sequence of the peptide, with modifications embedded.

        Returns:
        float: The total monoisotopic mass of the peptide, with modifications considered.

        Example:

        scales = Scales()
        print(scales.calculate_mass ("LVC[+57.02146]TALQW"))  # Modified peptide
        print(scales.calculate_mass_no_mods("LVC[+57.02146]TALQW"))  # Unmodified peptide

        Returns:
        972.473886
        915.4524260000001
        """
        # Regular expression to match modifications (e.g., "C[+57.02146]")
        pattern = re.compile(r'([A-Z])(\[\+?(-?\d+\.\d+)\])?')
        total_mass = 0.0

        for match in pattern.finditer(peptide):
            aa, _, mod_mass = match.groups()
            aa_mass = Scales.base_mass.get(aa, 0.0)
            total_mass += aa_mass + (float(mod_mass) if mod_mass else 0.0)
        
        return total_mass

    @staticmethod
    def calculate_mass_no_mods(peptide: str) -> float:
        """
        Calculate the total monoisotopic mass of an unmodified peptide sequence, ignoring any modifications.

        Parameters:
        peptide (str): The amino acid sequence of the peptide, with potential modifications embedded.

        Returns:
        float: The total monoisotopic mass of the peptide, with modifications ignored.


        Example:

        scales = Scales()
        print(scales.calculate_mass ("LVC[+57.02146]TALQW"))  # Modified peptide
        print(scales.calculate_mass_no_mods("LVC[+57.02146]TALQW"))  # Unmodified peptide

        Returns:
        972.473886
        915.4524260000001
        """
        # Remove modification notations from the peptide sequence
        clean_peptide = re.sub(r'\[.+?\]', '', peptide)
        # Calculate the mass using the base mass of the amino acids
        mass_list = [Scales.base_mass.get(aa, 0.0) for aa in clean_peptide]
        return sum(mass_list)
    
 
    @staticmethod
    def peptide_ipc(peptide: str, start_ph: float = 6.5, epsilon: float = 0.01) -> float:
        """
        Calculate the isoelectric point (pI) of a peptide.

        Parameters:
        peptide (str): The amino acid sequence of the peptide.
        start_ph (float): The starting pH value for pI calculation.
        epsilon (float): The precision for finding the pI value.

        Returns:
        float: The estimated isoelectric point of the peptide.

        # Example usage
        scales = Scales()
        print(scales.peptide_ipc("LVC[+57.02146]TALQW"))  # Modified peptide
    
        Returns:
        0.203125
        """
        IPC_score = {
            'Cterm': 2.383, 'pkD': 3.887, 'pkE': 4.317, 'pkC': 8.297, 'pkY': 10.071,
            'pkH': 6.018, 'Nterm': 9.564, 'pkK': 10.517, 'pkR': 12.503
        }

        peptide = re.sub(r'\[.+?\]', '', peptide)  # Strip modifications for IPC calculation

        aa_counts = {
            'D': peptide.count('D'),
            'E': peptide.count('E'),
            'C': peptide.count('C'),
            'Y': peptide.count('Y'),
            'H': peptide.count('H'),
            'K': peptide.count('K'),
            'R': peptide.count('R'),
        }

        nterm = peptide[0]
        cterm = peptide[-1]

        def charge_at_ph(ph_value, pk_value, is_positive):
            if is_positive:
                return 1.0 / (1.0 + 10 ** (ph_value - pk_value))
            else:
                return -1.0 / (1.0 + 10 ** (pk_value - ph_value))

        pH, pHprev, pHnext = start_ph, 0.0, 14.0

        while True:
            charge = (
                charge_at_ph(pH, IPC_score['Cterm'], cterm in ['D', 'E']) +
                charge_at_ph(pH, IPC_score['Nterm'], nterm in ['K', 'R', 'H']) +
                sum(charge_at_ph(pH, IPC_score[f'pk{aa}'], False) * aa_counts[aa] for aa in ['D', 'E', 'C', 'Y']) +
                sum(charge_at_ph(pH, IPC_score[f'pk{aa}'], True) * aa_counts[aa] for aa in ['H', 'K', 'R'])
            )

            if abs(charge) < epsilon:
                return pH

            if charge < 0.0:
                pHnext = pH
                pH -= (pH - pHprev) / 2.0
            else:
                pHprev = pH
                pH += (pHnext - pH) / 2.0

    @staticmethod
    def z_neutral_ph(peptide: str) -> float:
        """
        Calculate the net charge ('z') of a peptide at neutral pH.

        Basic amino acids (K, R) contribute a charge of +1, and acidic amino acids (D, E) contribute -1. 
        Other amino acids and modifications are not considered in this calculation.

        Parameters:
        peptide (str): The amino acid sequence of the peptide, potentially containing modifications.

        Returns:
        float: The net charge of the peptide at neutral pH.

        # Example usage:
        print(z_neutral_ph("LVC[+57.02146]TRLQW"))  # Modified peptide

        Expected Result:
        1.0 (assuming no basic or acidic amino acids in the sequence "LVC[+57.02146]TALQW")
        """
        z_dict = {'E': -1, 'D': -1, 'K': 1, 'R': 1}
        peptide = re.sub(r'\[.+?\]', '', peptide)  # Strip modifications, if any

        # Calculate net charge
        return sum(z_dict.get(aa, 0) for aa in peptide)

    @staticmethod
    def peptide_gravy(peptide: str, modifications: dict = None) -> float:
        """
        Calculate the Grand Average of Hydropathicity (GRAVY) of a peptide, considering modifications.

        The GRAVY value is calculated as the sum of hydropathy values of all the amino acids and modifications 
        in the peptide divided by the number of residues in the sequence.

        Parameters:
        peptide (str): The amino acid sequence of the peptide, potentially containing modifications.
        modifications (dict, optional): A dictionary of modified amino acids and their GRAVY scores.

        Returns:
        float: The GRAVY score of the peptide.

        # Example usage:
        print(peptide_gravy("ALWKTMLKY"))  # Unmodified peptide
        print(peptide_gravy("ALWKC[+57.6857]TMLKY", {"C_57_6857": 3.2}))  # Modified peptide

        Expected Results:
        0.06  #Unmodified
        0.38  #Modified
        """
        base_hydro = {
        "A": 1.800, "R": -4.500, "N": -3.500, "D": -3.500, "C": 2.500, "Q": -3.500,
            "E": -3.500, "G": -0.400, "H": -3.200, "I": 4.500, "L": 3.800, "K": -3.900,
            "M": 1.900, "F": 2.800, "P": -1.600, "S": -0.800, "T": -0.700, "W": -0.900,
            "Y": -1.300, "V": 4.200,
        }

        if modifications is None:
            modifications = {}

        # Function to parse and get the GRAVY score of each amino acid
        def get_hydro_score(aa):
            if '[' in aa and ']' in aa:
                mod_key = aa[0] + '_' + aa[aa.find('[') + 1: aa.find(']')].replace('+', '').replace('.', '_')
                return modifications.get(mod_key, 0.0)
            return base_hydro.get(aa, 0.0)

        # Split the peptide into amino acids and modifications
        peptide_components = re.findall(r'[A-Z]\[\+\d+\.\d+\]|[A-Z]', peptide)

        if len(peptide_components) == 0:
            return 0.0  # Return 0.0 to avoid division by zero

        hydro_scores = [get_hydro_score(aa) for aa in peptide_components]
        return sum(hydro_scores) / len(peptide_components)
    
    @staticmethod
    def peptide_inspector(peptide: str, percentage: bool = False) -> dict:
        """
        Analyze the composition of a peptide and return the count or percentage of each amino acid.

        Parameters:
        peptide (str): The amino acid sequence of the peptide.
        percentage (bool): If True, return the percentage composition of each amino acid; otherwise, return the raw count.

        Returns:
        dict: A dictionary with amino acids as keys and their counts or percentages as values.

        # Example usage:
        print(peptide_inspector("ALWKTMLKY"))  # Raw count
        print(peptide_inspector("ALWKTMLKY", percentage=True))  # Percentage

        Expected Results:
        {"A": 1, "L": 2, "W": 1, "K": 2, "T": 1, "M": 1, "Y": 1}  # Raw count
        {"A": 0.11, "L": 0.22, "W": 0.11, "K": 0.22, "T": 0.11, "M": 0.11, "Y": 0.11}  # Percentage
        """
        aa_count = dict.fromkeys(set(peptide), 0)
        for aa in peptide:
            aa_count[aa] += 1

        if percentage:
            total_length = len(peptide)
            for aa in aa_count:
                aa_count[aa] = round(aa_count[aa] / total_length, 2)

        return aa_count


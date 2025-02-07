import random
import re


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
    def calculate_peptide_charge(peptide, pH):
        # pKa values for the N-terminus, C-terminus, and side chains of ionizable amino acids
        pKa = {
            'N_term': 9.69,
            'C_term': 2.34,
            'K': 10.4,
            'R': 12.5,
            'H': 6.0,
            'D': 3.9,
            'E': 4.1,
            'C': 8.3,
            'Y': 10.1,
        }

        # Charge contributions by amino acid at the given pH
        charge = {
            'N_term': 1 / (1 + pow(10, pH - pKa['N_term'])),
            'C_term': -1 / (1 + pow(10, pKa['C_term'] - pH)),
            'K': 1 / (1 + pow(10, pH - pKa['K'])),
            'R': 1 / (1 + pow(10, pH - pKa['R'])),
            'H': 1 / (1 + pow(10, pH - pKa['H'])),
            'D': -1 / (1 + pow(10, pKa['D'] - pH)),
            'E': -1 / (1 + pow(10, pKa['E'] - pH)),
            'C': -1 / (1 + pow(10, pKa['C'] - pH)),
            'Y': -1 / (1 + pow(10, pKa['Y'] - pH)),
        }

        # Calculate the net charge
        net_charge = charge['N_term'] + charge['C_term']
        for aa in peptide:
            if aa in charge:
                net_charge += charge[aa]

        return round(net_charge)


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


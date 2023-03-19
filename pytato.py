#Pytato 2023 

import argparse
import csv
import os
import pandas as pd
from pyteomics import fasta, mass, parser
import subprocess
from support import *

#````````````````````````````````````````````````````````````````````````````````````````````````````````


##Stage 1: Convert .RAW files to .mzML files using Proteowizard/MSConvert
##

def convert_raw_to_mzml(input_folder, output_folder, msconvert_path):
    raw_files = [f for f in os.listdir(input_folder) if f.lower().endswith('.raw')]
    os.makedirs(output_folder,exist_ok=True)

    if not raw_files:
        print("No .RAW files found in the input folder.")
        return

    for raw_file in raw_files:
        input_file = os.path.join(input_folder, raw_file)
        print(f"Converting {input_file} to mzML...")
        cmd = f"{msconvert_path} {input_file} -o {output_folder} --mzML"
        subprocess.run(cmd, shell=True, check=True)

# # Example usage
# input_folder = "/path/to/input_folder"
# output_folder = "/path/to/output_folder"
# msconvert_path = "/path/to/msconvert.exe"


#````````````````````````````````````````````````````````````````````````````````````````````````````````

## STAGE 2: SPECTRAL LIBRARY BUILD WITH DIA-Umpire
## 

def generate_spectral_library(input_folder, output_folder, dia_umpire_jar_path):
    mzml_files = [f for f in os.listdir(input_folder) if f.lower().endswith('.mzml')]
    os.makedirs(output_folder,exist_ok=True)
    if not mzml_files:
        print("No .mzML files found in the input folder.")
        return

    # Set the parameters for DIA-Umpire
    params = {
        'AdjustFragIntensity': 'true',
        'MS1PPM': '10',
        'MS2PPM': '20',
        'MS2SN': '1.5',
        'MS1SN': '1.5',
        'Thread': '4',
        'Ram': '20000',
        'Quant': 'false'
    }

    # Convert the parameters dictionary to a string format
    params_str = ' '.join([f"-{k} {v}" for k, v in params.items()])

    for mzml_file in mzml_files:
        input_file = os.path.join(input_folder, mzml_file)
        print(f"Generating spectral library for {input_file}...")
        
        # Build the command to run DIA-Umpire
        cmd = f"java -jar {dia_umpire_jar_path} {input_file} -{output_folder} {params_str}"
    
        # Execute the command
        subprocess.run(cmd, shell=True, check=True)

# Example usage
# input_file = "/path/to/your/mzML_file.mzML"
# output_folder = "/path/to/output_folder"
# dia_umpire_jar_path = "/path/to/DIA_Umpire_SE.jar"


#````````````````````````````````````````````````````````````````````````````````````````````````````````


##Stage 3: DIA-NN Search
##

def run_dia_nn(input_folder, library_file, output_folder, dia_nn_exe_path):
    mzml_files = [f for f in os.listdir(input_folder) if f.lower().endswith('.mzml')]

    if not mzml_files:
        print("No .mzML files found in the input folder.")
        return

    for mzml_file in mzml_files:
        input_file = os.path.join(input_folder, mzml_file)
        output_file = os.path.join(output_folder, f"{os.path.splitext(mzml_file)[0]}_output.tsv")
        print(f"Processing {input_file} with DIA-NN...")

        # Build the command to run DIA-NN
        cmd = f"{dia_nn_exe_path} --f {input_file} --lib {library_file} --out {output_file}"

        # Execute the command
        subprocess.run(cmd, shell=True, check=True)

# # Example usage
# input_folder = "/path/to/input_folder"
# library_file = "/path/to/library_file"
# output_folder = "/path/to/output_folder"
# dia_nn_exe_path = "/path/to/dia_nn.exe"

## Step 4: Produce Theoretical Spectra for 2nd Search using High-Confidence Hits from 1st Search
#
def parse_dia_nn_results(result_file, confidence_threshold):
    results_df = pd.read_csv(result_file, sep='\t')
    high_conf_proteins = results_df[results_df['Confidence'] >= confidence_threshold]['Protein.Id'].unique()
    return high_conf_proteins

result_file = "/path/to/dia_nn_result.tsv"
confidence_threshold = 0.99

high_conf_proteins = parse_dia_nn_results(result_file, confidence_threshold)



def generate_digested_peptides(fasta_file, high_conf_proteins, cleavage_rule):
    with fasta.read(fasta_file) as f:
        for entry in f:
            protein_id = entry.description
            protein_seq = entry.sequence
            if protein_id in high_conf_proteins:
                peptides = parser.cleave(protein_seq, cleavage_rule)
                for peptide in peptides:
                    yield peptide

fasta_file = "/path/to/protein_database.fasta"
thermolysin_cleavage_rule = '[ALIVGFYW]'

digested_peptides = list(generate_digested_peptides(fasta_file, high_conf_proteins, thermolysin_cleavage_rule))


def generate_theoretical_spectra(peptides, charge_range=(1, 3), output_file="theoretical_spectra.tsv"):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Peptide', 'Charge', 'IonType', 'IonNumber', 'mz'])

        for peptide in peptides:
            for charge in range(charge_range[0], charge_range[1] + 1):
                spectrum = mass.generate_spectrum(peptide, charge=charge, ion_types=('b', 'y'))
                for ion_type, ion_series in spectrum.items():
                    for ion_number, mz in enumerate(ion_series, start=1):
                        writer.writerow([peptide, charge, ion_type, ion_number, mz])






#````````````````````````````````````````````````````````````````````````````````````````````````````````

def main():
    parser = argparse.ArgumentParser(description='Setup Pytato Enviroment')
    parser.add_argument('--pytato-folder', required=True, help='Path to the Pytato folder.')
    parser.add_argument('--input_folder',default=os.listdir(),help='Folder containing .RAW files')
    parser.add_argument('--mzml_folder',default=os.listdir(),help='Folder containing mzml files')
    parser.add_argument('--output_folder',default=os.listdir(),help='Output Folder')
    parser.add_argument('--dia-umpire-url', default='https://github.com/diaumpire/DIA-Umpire/releases/download/v2.2/DIA_Umpire_SE.jar', help='Direct download URL for DIA_Umpire_SE.jar from GitHub (default: %(default)s).')
    parser.add_argument('--dia_umpire_jar_path', default="No file path specified",help="Path to dia_umpire_jar_path")
    parser.add_argument('--dia_nn_exe_path', default="No file path specified",help="Path to dia_nn .exe")
    parser.add_argument('--msconvert_path', default="No file path specified",help="Path to msconvert .exe")
    parser.add_argument('--pull_msconvert', default='N')
    parser.add_argument('--pull_diann', default='N')
    parser.add_argument('--pull_umpire', default='N')
    parser.add_argument('--bake', default='OFF')

    args = parser.parse_args()
    if args.pull_msconvert == "Y":
        download_msconvert(args.pytato_folder)
    if args.pull_umpire == "Y":
        download_dia_umpire(args.pytato_folder, args.dia_umpire_url)
    if args.pull_diann == "Y":
        download_dia_nn(args.pytato_folder,args.dia_nn_url)
    
    if args.bake == "ON":
        #Preheat
        convert_raw_to_mzml(args.input_folder, args.mzml_folder, args.msconvert_path) ## RAW to mzML
        generate_spectral_library(args.output_folder, args.output_folder, args.dia_umpire_jar_path)
        #First Bake/Search (Forward)
        run_dia_nn(args.input_folder, args.lib_file1, args.output_folder, args.dia_nn_exe_path)
        #Second Bake/Search (Forward)
        generate_theoretical_spectra(digested_peptides, output_file="theoretical_spectra.tsv")

        #First Search (Reverse)
        #Second Search (Reverse)
        #Concatenate/Summarize Data




if __name__ == '__main__':
    main()
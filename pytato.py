#Pytato 2023 

import argparse
import csv
import os
import pandas as pd
from pyteomics import fasta, mass, parser
import subprocess
from support import *


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


def generate_spectral_library(input_folder, output_folder, dia_umpire_jar_path):
    mzml_files = [f for f in os.listdir(input_folder) if f.lower().endswith('.mzml')]
    os.makedirs(output_folder, exist_ok=True)
    if not mzml_files:
        print("No .mzML files found in the input folder.")
        return []

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

    mgf_files = []

    for mzml_file in mzml_files:
        input_file = os.path.join(input_folder, mzml_file)
        print(f"Generating spectral library for {input_file}...")

        # Build the command to run DIA-Umpire
        cmd = f"java -jar {dia_umpire_jar_path} {input_file} -{output_folder} {params_str}"

        # Execute the command
        subprocess.run(cmd, shell=True, check=True)

        mgf_file = os.path.join(output_folder, f"{os.path.splitext(mzml_file)[0]}_Q1.mgf")
        mgf_files.append(mgf_file)

    return mgf_files



def fasta_to_proteins(fasta_file):
    with open(fasta_file, 'r') as file_handle:
        content = file_handle.read()

    proteins = []
    fasta_entries = content.split('>')[1:]  # Remove the first empty entry
    for entry in fasta_entries:
        lines = entry.strip().split('\n')
        header, sequence = lines[0], ''.join(lines[1:])
        proteins.append(sequence)

    return proteins


def run_dia_nn(input_folder, library_files, output_folder, dia_nn_exe_path):
    mzml_files = [f for f in os.listdir(input_folder) if f.lower().endswith('.mzml')]

    if not mzml_files:
        print("No .mzML files found in the input folder.")
        return []

    report_files = []

    for mzml_file in mzml_files:
        input_file = os.path.join(input_folder, mzml_file)
        output_file_prefix = os.path.join(output_folder, os.path.splitext(mzml_file)[0])
        report_file = f"{output_file_prefix}.report.tsv"
        report_files.append(report_file)
        print(f"Processing {input_file} with DIA-NN...")

        library_str = ' '.join([f"--lib {lib}" for lib in library_files])

        # Build the command to run DIA-NN
        cmd = f"{dia_nn_exe_path} --f {input_file} {library_str} --out {output_file_prefix}"

        # Execute the command
        subprocess.run(cmd, shell=True, check=True)

    return report_files

def get_high_confidence_proteins(report_tsv, fdr_threshold=0.01):
    df = pd.read_csv(report_tsv, sep='\t')
    high_conf_proteins = df[df['Protein.Q.Value'] <= fdr_threshold]['Protein.Ids']
    return set(high_conf_proteins)


def generate_digested_peptides(fasta_file, proteins, cleavage_rule):
    fasta_proteins = fasta_to_proteins(fasta_file)
    for protein_id, protein_seq in fasta_proteins.items():
        if protein_id in proteins:
            peptides = parser.cleave(protein_seq, cleavage_rule)
            for peptide in peptides:
                yield peptide

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
    
    return output_file






#````````````````````````````````````````````````````````````````````````````````````````````````````````

def main():
    parser = argparse.ArgumentParser(description='Setup Pytato Enviroment')
    #args to prepare search engine
    parser.add_argument('--pytato-folder', required=True, help='Path to the Pytato folder.')
    parser.add_argument('--input_folder1',default=os.listdir(),help='Folder containing .RAW files for Enzyme1')
    parser.add_argument('--input_folder2',default=os.listdir(),help='Folder containing .RAW files for Enzyme2')
    parser.add_argument('--mzml_folder',default=os.listdir(),help='Folder containing mzml files')
    parser.add_argument('--output_folder',default=os.listdir(),help='Output Folder')
    parser.add_argument('--fasta_file_path',default="No FASTA File Provided",help='Download FASTA file and provide path')
    parser.add_argument('--dia-umpire-url', default='https://github.com/diaumpire/DIA-Umpire/releases/download/v2.2/DIA_Umpire_SE.jar', help='Direct download URL for DIA_Umpire_SE.jar from GitHub (default: %(default)s).')
    parser.add_argument('--dia_umpire_jar_path', default="No file path specified",help="Path to dia_umpire_jar_path")
    parser.add_argument('--dia_nn_exe_path', default="No file path specified",help="Path to dia_nn .exe")
    parser.add_argument('--msconvert_path', default="No file path specified",help="Path to msconvert .exe")
    parser.add_argument('--pull_msconvert', default='N')
    parser.add_argument('--pull_diann', default='N')
    parser.add_argument('--pull_umpire', default='N')
    #args to perform search
    parser.add_argument('--bake', default='OFF', help="Turn ON if to perform search")
    parser.add_argument('--enzyme1_name', default='Enzyme1',help='User-defined Name for Enzyme1')
    parser.add_argument('--enzyme1_rule', default='trypsin',help='Cleavage Rule for Enzyme1')
    parser.add_argument('--enzyme1_mc', default='2',help='Missed Clevage Number for Enzyme1')
    parser.add_argument('--enzyme2_name', default='Enzyme2',help='User-defined Name for Enzyme2')    
    parser.add_argument('--enzyme2_rule', default='thermolysin',help='Cleavage Rule for Enzyme2')
    parser.add_argument('--enzyme2_mc', default='2',help='Missed Clevage Number for Enzyme2')
    parser.add_argument('--confidence_lvl_1', default="0.90",help="Confidence level for Forward Search (First Bake), Selects Proteins with high-confidence from 1st DIA-NN Seach (Enzyme1)")
    parser.add_argument('--confidence_lvl_2', default="0.90",help="Confidence level for Reverse Search (Second Bake), Selects Proteins with high-confidence from 2nd DIA-NN Seach (Enzyme2)")
    
    args = parser.parse_args()
    ## ENVIROMENT SETUP
    if args.pull_msconvert == "Y":
        download_msconvert(args.pytato_folder)
    if args.pull_umpire == "Y":
        download_dia_umpire(args.pytato_folder, args.dia_umpire_url)
    if args.pull_diann == "Y":
        download_dia_nn(args.pytato_folder,args.dia_nn_url)
    
    if args.bake == "ON": ## EXECUTE SEARCH
        ##Preheat
        convert_raw_to_mzml(args.input_folder1, args.mzml_folder, args.msconvert_path) ## RAW to mzML
        mgf1=generate_spectral_library(args.output_folder, args.output_folder, args.dia_umpire_jar_path)

        ##First Bake/Search (Forward)
        total_proteins=fasta_to_proteins(args.fasta_file_path)
        peptide_generator = generate_digested_peptides(args.fasta_file_path, total_proteins, args.enzyme1_rule)
        spectra1=generate_theoretical_spectra(peptide_generator)        
        results=run_dia_nn(mgf1, spectra1, args.output_folder, args.dia_nn_exe_path)
        confidence_threshold1 = int(args.confidence_lvl_1)
        high_conf_proteins=[]
        for result in results:
            highcon_pro = get_high_confidence_proteins(result, confidence_threshold1)
            high_conf_proteins.append(highcon_pro)
        high_conf_proteins=list(set(high_conf_proteins))
        ##Second Bake/Search (Forward)
        peptide_generator = generate_digested_peptides(args.fasta_file_path, high_conf_proteins, args.enzyme2_rule)
        spectra2=generate_theoretical_spectra(peptide_generator, output_file=f"{args.enzyme2_name}_spectra.tsv")
        results=run_dia_nn(args.input_folder, spectra2, args.output_folder, args.dia_nn_exe_path)


        #First Search (Reverse)
        #Second Search (Reverse)
        #Concatenate/Summarize Data



if __name__ == '__main__':
    main()
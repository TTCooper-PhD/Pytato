#Pytato 2023 

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import os
import pandas as pd
from pathlib import Path
from pyteomics import parser,fasta,mass, mgf,auxiliary
import shutil
import subprocess
from support import *


def convert_raw_to_mzml(input_folder, msconvert_path, output_subfolder="mzML"):
    """
    Converts .RAW files in the input folder to .mzML files and saves them in a subdirectory.
    
    Parameters:
    input_folder (str): Path to the folder containing .RAW files.
    msconvert_path (str): Path to the msconvert executable.
    output_subfolder (str, optional): Name of the subdirectory to save the converted .mzML files. Defaults to "mzML".
    
    Returns:
    str: Path to the output folder where the converted .mzML files are saved.
    """
    raw_files = [f for f in os.listdir(input_folder) if f.lower().endswith('.raw')]

    if not raw_files:
        print("No .RAW files found in the input folder.")
        return None

    output_folder = os.path.join(input_folder, output_subfolder)
    os.makedirs(output_folder, exist_ok=True)

    for raw_file in raw_files:
        input_file = os.path.join(input_folder, raw_file)
        print(f"Converting {input_file} to mzML...")
        cmd = f"\"{msconvert_path}\" {input_file} -o {output_folder} --mzML --filter \"peakPicking vendor msLevel=1\" --filter \"zeroSamples removeExtra\""
        subprocess.run(cmd, shell=True, check=True)

    return output_folder


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

def generate_spectral_library(dia_nn_exe_path, fasta_file):
    """
    Generates a spectral library from a given FASTA file using DIA-NN.

    Args:
        dia_nn_exe_path (str): Path to the DIA-NN executable.
        fasta_file (str): Path to the FASTA file.

    Returns:
        str: Path to the generated spectral library file ('spec-lib-predicted.dlib').
    """
    if not os.path.exists(fasta_file):
        print(f"FASTA file '{fasta_file}' not found.")
        return ""

    output_file = "spec-lib-predicted.dlib"

    cmd = f"{dia_nn_exe_path} \
    --lib \"\" \
    --threads 30 \
    --verbose 3 \
    --out \"report.tsv\" \
    --qvalue 0.01 \
    --matrices \
    --out-lib \"report-lib.tsv\" \
    --gen-spec-lib \
    --predictor \
    --fasta \"{fasta_file}\" \
    --fasta-search"

    # Execute the command
    subprocess.run(cmd, shell=True, check=True)

    if os.path.exists(output_file):
        return list(os.path.abspath(output_file))
    else:
        print("Spectral library file not found.")
        return ""
    
def run_dia_nn(dia_nn_exe_path, library_files=list(), fasta_files=list(), input_folder=str(), output_folder=str(),report_file_name=str(), sn_ratio=1.0,qval=0.01, threads=30,missed_cleavages=1,
               cut="K*,R*",min_frag_mz=200,max_frag_mz=1800,min_pre_mz=300,max_pre_mz=1200, min_pep_len=7,max_pep_len=30,
               min_pre_z=1,max_pre_z=4):
    """
    Runs DIA-NN on a set of .mzML files from a specified input folder using specified parameters.

    Args:
        dia_nn_exe_path (str): Path to the DIA-NN executable.
        library_files (list): List of spectral library files (.speclib).
        fasta_files (list): List of FASTA files (.fasta).
        input_folder (str): Path to the input folder containing .mzML files.
        output_folder (str): Path to the output folder where the results will be saved.
        report_file_name (str): Name of the output report file (without the extension).
        sn_ratio (float): Signal-to-noise ratio.
        qval (float): Target q-value.
        threads (int): Number of threads to use for processing.
        missed_cleavages (int): Maximum number of missed cleavages allowed.
        cut (str): Protease cleavage rule.
        min_frag_mz (int): Minimum fragment m/z.
        max_frag_mz (int): Maximum fragment m/z.
        min_pre_mz (int): Minimum precursor m/z.
        max_pre_mz (int): Maximum precursor m/z.
        min_pep_len (int): Minimum peptide length.
        max_pep_len (int): Maximum peptide length.
        min_pre_z (int): Minimum precursor charge.
        max_pre_z (int): Maximum precursor charge.

    Returns:
        str: Path to the generated report file.
    """
    mzml_files = [f"/{f}" for f in os.listdir(input_folder) if f.lower().endswith('.mzml')]
    os.makedirs(output_folder,exist_ok=True)
    if not mzml_files:
        print("No .mzML files found in the input folder.")
        return []

    file_str = ' '.join([f"--lib {fil}" for fil in mzml_files])
    library_str = ' '.join([f"--lib {lib}" for lib in library_files])
    fasta_str = ' '.join([f"--fasta {fasta}" for fasta in fasta_files])
    report_file=f"'{output_folder}/{report_file_name}.tsv'"

    cmd = f"{dia_nn_exe_path} \
    {file_str} \
    {library_str}\
    {threads} \
    --verbose 3 \
    --out {report_file} \
    --qvalue {qval} \
    --matrices \
    --out-lib f'{report_file_name}-lib.tsv\" \
    --gen-spec-lib \
    --predictor \
    {fasta_str}\
    --fasta-search \
    --min-fr-mz {min_frag_mz} \
    --max-fr-mz {max_frag_mz} \
    --met-excision --cut {cut} \
    --missed-cleavages {missed_cleavages} \
    --min-pep-len {min_pep_len} --max-pep-len {max_pep_len} \
    --min-pr-mz {min_pre_mz} --max-pr-mz {max_pre_mz} \
    --min-pr-charge {min_pre_z} --max-pr-charge {max_pre_z} \
    --unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n --monitor-mod UniMod:1"

    # Execute the command
    subprocess.run(cmd, shell=True, check=True)

    return report_file


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

def concatenate_results(results_1, enzyme1_name, results_2, enzyme2_name):
    # Read search results as dataframes
    df_1 = pd.read_csv(results_1)
    df_2 = pd.read_csv(results_2)

    # Create a new column 'Detected_in' and set values to enzyme name
    df_1['Detected_in'] = f'{enzyme1_name}'
    df_2['Detected_in'] = f'{enzyme2_name}'

    # Concatenate dataframes
    combined_results = pd.concat([df_1, df_2], axis=0, ignore_index=True)

    # Create a new column 'Detected_in_both' and set its default value to False
    combined_results['Detected_in_both'] = False

    # Group the combined dataframe by protein and iterate through groups
    grouped = combined_results.groupby('Protein')
    for protein, group in grouped:
        # Check if protein is detected in both enzymes
        if group['Detected_in'].nunique() == 2:
            # Update 'Detected_in_both' column for the protein
            combined_results.loc[group.index, 'Detected_in_both'] = True

    return combined_results

def run_search(args,direction):
    raw_files = [f for f in os.listdir(args.input_folder) if f.lower().endswith('.raw')]
    
    if raw_files:
        convert_raw_to_mzml(args.input_folder, args.mzml_folder, args.msconvert_path)  # RAW to mzML


    if direction == "forward":
        ##First Bake/Search (Forward)
        total_proteins=fasta_to_proteins(args.fasta_file_path) #pull down proteins from FASTA file
        peptide_generator = generate_digested_peptides(args.fasta_file_path, total_proteins, args.enzyme1_rule) #generate peptides from total_proteins using enzyme1
        spectra1=generate_theoretical_spectra(peptide_generator) #generate theoretical spectrum with peptides        
        results=run_dia_nn(mgf1, spectra1, args.output_folder, args.dia_nn_exe_path) #run DIA-NN using mgf files (Group1)
        confidence_threshold1 = int(args.confidence_lvl_1) #Determin Q-value threshold
        high_conf_proteins=[] #Collect Proteins Identified with High Confidence (e.g. FDR=0.01)
        for result in results:
            #Collect Proteins Identified with High Confidence from Each Sample Analyzed
            highcon_pro = get_high_confidence_proteins(result, confidence_threshold1)
            high_conf_proteins.append(highcon_pro)
        high_conf_proteins=list(set(high_conf_proteins))
        ##Second Bake/Search (Forward)
        peptide_generator = generate_digested_peptides(args.fasta_file_path, high_conf_proteins, args.enzyme2_rule) #Use list of high-confidence proteins from first search to refine theoretical peptides for 2nd search
        spectra2=generate_theoretical_spectra(peptide_generator, output_file=f"{args.enzyme2_name}_spectra.tsv") #Generare theoretical spectra library
        output_fwd=f"{args.output}/Output_{args.enzyme2_name}"
        results=run_dia_nn(mgf2, spectra2, args.output_folder, args.dia_nn_exe_path) #Run DIA-NN

    elif direction == "reverse":
        #First Search (Reverse)
        peptide_generator = generate_digested_peptides(args.fasta_file_path, total_proteins, args.enzyme2_rule) #generate peptides from total_proteins using enzyme1
        spectra1=generate_theoretical_spectra(peptide_generator) #generate theoretical spectrum with peptides        
        results=run_dia_nn(mgf2, spectra1, args.output_folder, args.dia_nn_exe_path) #run DIA-NN using mgf files (Group1)
        confidence_threshold1 = int(args.confidence_lvl_1) #Determin Q-value threshold
        high_conf_proteins=[] #Collect Proteins Identified with High Confidence (e.g. FDR=0.01)
        for result in results:
            #Collect Proteins Identified with High Confidence from Each Sample Analyzed
            highcon_pro = get_high_confidence_proteins(result, confidence_threshold1)
            high_conf_proteins.append(highcon_pro)

        #Second Search (Reverse)
        peptide_generator = generate_digested_peptides(args.fasta_file_path, high_conf_proteins, args.enzyme1_rule) #Use list of high-confidence proteins from first search to refine theoretical peptides for 2nd search
        spectra2=generate_theoretical_spectra(peptide_generator, output_file=f"{args.enzyme1_name}_spectra.tsv") #Generare theoretical spectra library
        output_rev=f"{args.output}/Output_{args.enzyme1_name}"
        results=run_dia_nn(mgf1, spectra2, output_rev, args.dia_nn_exe_path) #Run DIA-NN

    else:
        raise ValueError("Invalid search direction.")


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
    parser.add_argument('--dia_nn_exe_path', default="No file path specified",help="Path to dia_nn .exe")
    parser.add_argument('--msconvert_path', default="No file path specified",help="Path to msconvert .exe")
    parser.add_argument('--pull_msconvert', default='N')
    parser.add_argument('--pull_diann', default='N')
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
    parser.add_argument('--sn_ratio', default="1", help='Signal-to-noise ratio threshold for DIA-NN search (default: %(default)s).')   
    args = parser.parse_args()
    ## ENVIROMENT SETUP


if __name__ == '__main__':
    main()
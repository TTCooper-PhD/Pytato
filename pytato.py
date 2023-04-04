#Pytato 2023 

import argparse
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
import csv
import os
import pandas as pd
from pathlib import Path
from pyteomics import parser,fasta,mass, mgf,auxiliary
import re
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
        cmd = f"\"{msconvert_path}\" {input_file} -o {output_folder} --mzML --filter \"peakPicking vendor msLevel=1\" --filter \"zeroSamples removeExtra 1\" --filter \"demultiplex optimization=overlap massError=10ppm\""
        subprocess.run(cmd, shell=True, check=True)

    return output_folder



def filter_fasta(input_fasta, output_fasta, uniprot_ids):
    """
    Filters a FASTA file to include only sequences with the specified UniProt IDs.

    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to the output filtered FASTA file.
        uniprot_ids (list): List of UniProt IDs to keep in the output FASTA file.

    Returns:
        str: Path to Output FASTA
    """
    uniprot_id_pattern = re.compile(r"sp\|(\w+)\|")

    with open(input_fasta, "r") as in_file, open(output_fasta, "w") as out_file:
        fasta_sequences = SeqIO.parse(in_file, "fasta")
        filtered_sequences = []

        for seq in fasta_sequences:
            match = uniprot_id_pattern.search(seq.description)
            if match:
                uniprot_id = match.group(1)
                if uniprot_id in uniprot_ids:
                    filtered_sequences.append(seq)

        SeqIO.write(filtered_sequences, out_file, "fasta")

    return output_fasta


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
    --fasta-search \
    --smart-profiling\
    --peak-center \
    --no-ifs-removal"

    # Execute the command
    subprocess.run(cmd, shell=True, check=True)

    if os.path.exists(output_file):
        return os.path.abspath(output_file)
    else:
        print("Spectral library file not found.")
        return ""
    


    
def run_dia_nn(dia_nn_exe_path, library_files, fasta_files, input_folder, output_folder,report_file_name="report", qval=0.01,threads=30,missed_cleavages=1,
               cut="K*,R*",min_frag_mz=200,max_frag_mz=1800,min_pre_mz=300,max_pre_mz=1200, min_pep_len=7,max_pep_len=30,ms2_acc=20,ms2_acc_cal=20,ms1_acc=20,
               min_pre_z=1,max_pre_z=4,fasta_search=False,profiling="smart",MBR=True,fasta_speclib_annotation=False,frag_restrict_quant=True,heuristic_search=True):

    mzml_files = [f"{f}" for f in os.listdir(input_folder) if f.lower().endswith('.mzml')]
    os.makedirs(output_folder,exist_ok=True)
    if not mzml_files:
        print("No .mzML files found in the input folder.")
        return []
    file_str = ' '.join([f"--f {fil}" for fil in mzml_files])

    if isinstance(library_files, str):
        library_files = [library_files]
    library_str = ' '.join([f"--lib {lib}" for lib in library_files])

    if isinstance(fasta_files, str):
        fasta_files = [fasta_files]
    fasta_str = ' '.join([f"--fasta {fasta}" for fasta in fasta_files])

    report_file=f"{output_folder}/{report_file_name}.tsv"

    if fasta_search==True:
        fasta_search_out="--fasta-search"
    else:
        fasta_search_out=""
    
    if MBR==True:
        match_between_runs="--reanalyze"
    else:
        match_between_runs=""    

    if fasta_speclib_annotation==True:
        reannotate="--reannotate"
    else:
        reannotate=""   

    if frag_restrict_quant==True:
        fr_r_quant="--gen-fr-restriction"
    else:
        fr_r_quant=""  
           
    
    if heuristic_search==True:
        relaxed_prot_inf="--relaxed-prot-inf"
    else:
        relaxed_prot_inf=""  

    if profiling.lower()=="smart":
        profiling_out="--smart-profiling"
    elif profiling.lower()=="rt_profiling":
        profiling_out=="--rt-profiling"
    else:
        profiling_out="--smart-profiling"

    cmd = f"{dia_nn_exe_path} \
    {file_str} \
    {library_str}\
    --threads {threads} \
    --verbose 4 \
    --out {report_file} \
    --qvalue {qval} \
    --matrices \
    --out-lib {report_file_name}-lib.tsv \
    --gen-spec-lib \
    --predictor \
    {fasta_str}\
    {fasta_search_out} \
    --min-fr-mz {min_frag_mz} \
    --max-fr-mz {max_frag_mz} \
    --met-excision \
    --cut {cut} \
    --missed-cleavages {missed_cleavages} \
    --min-pep-len {min_pep_len} --max-pep-len {max_pep_len} \
    --min-pr-mz {min_pre_mz} --max-pr-mz {max_pre_mz} \
    --min-pr-charge {min_pre_z} --max-pr-charge {max_pre_z} \
    --mass-acc {ms2_acc} \
    --mass-acc-cal {ms2_acc_cal}\
    --mass-acc-ms1 {ms1_acc}\
    {match_between_runs} \
    {reannotate}\
    {fr_r_quant} \
    --peak-center \
    {profiling_out} \
    {heuristic_search} \
    --no-ifs-removal\
    --unimod4 \
    --var-mods 1 \
    --var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n \
    --monitor-mod UniMod:1"

    # Execute the command
    subprocess.run(cmd, shell=True, check=True)

    return report_file


def get_high_confidence_proteins(report_tsv, fdr_threshold=0.01):
    df = pd.read_csv(report_tsv, sep='\t')
    high_conf_proteins = df[df['Protein.Q.Value'] <= fdr_threshold]['Protein.Ids']
    return set(high_conf_proteins)



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

    elif direction == "reverse":

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
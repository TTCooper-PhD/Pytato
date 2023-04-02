from pytato import *

#Step 1:
#Generate mzML files from .RAW files
raw_files="Target_Folder"
os.listdir("C:\Program Files\ProteoWizard\ProteoWizard 3.0.23")
msconvert_path="C:\Program Files\ProteoWizard\ProteoWizard 3.0.23\msconvert"
mzml_folder=convert_raw_to_mzml(raw_files,msconvert_path,output_subfolder="mzML")
mzml_folder="Target_Folder\mzML"
#Step 2:
#Generate Spectral Library Using mzML files Generated in Step 1
dia_umpire_path="DIA_Umpire_SE.jar"
mgf_files=generate_spectral_library(mzml_folder,"mgf_files",dia_umpire_path)
diann_input_folder="diann_input"
files_tsv_4_diann=mgf_to_tsv(mgf_files,diann_input_folder)

#Step 3: Generate List of Candidate Proteins
fasta_file="H:\\Tyler_Cooper\\Uniprot_Human_Curated_July_2020.fasta"
total_proteins=fasta_to_proteins(fasta_file)
print(len(total_proteins))

#Step 4: Run DIA-NN (First Bake - Forward)
diann_path="C:\\DIA-NN\\1.8.1\\diaNN.exe"
target_files = [os.path.join(f"{diann_input_folder}", f) for f in os.listdir(f"{diann_input_folder}") if f.lower().endswith(".tsv")]
report_files = run_dia_nn(target_files, mzml_folder, "FirstBake_Forward", diann_path, 1.0)


cmd = f"{diann_path} \
--f \"Target_Folder\File1_-35.mzML\" \
--lib \"\" \
--threads 30 \
--verbose 3 \
--out \"report.tsv\" \
--qvalue 0.01 \
--matrices \
--out-lib \"report-lib.tsv\" \
--gen-spec-lib \
--predictor \
--fasta \"FASTA_HUMAN_Apr2023.fasta\" \
--fasta-search \
--min-fr-mz 200 \
--max-fr-mz 1800 \
--met-excision --cut K*,R* \
--missed-cleavages 1 \
--min-pep-len 7 --max-pep-len 30 \
--min-pr-mz 300 --max-pr-mz 1800 \
--min-pr-charge 1 --max-pr-charge 4 \
--unimod4 --var-mods 1 --var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n --monitor-mod UniMod:1"

# Execute the command
subprocess.run(cmd, shell=True, check=True)
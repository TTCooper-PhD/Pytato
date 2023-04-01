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

#Step 3: Generate List of Candidate Proteins
fasta_file="H:\\Tyler_Cooper\\Uniprot_Human_Curated_July_2020.fasta"
total_proteins=fasta_to_proteins(fasta_file)
print(len(total_proteins))

#Step 4: Run DIA-NN (First Bake - Forward)
diann_path="C:\\DIA-NN\\1.8.1\\diaNN.exe"
mgf_files = [os.path.join("mgf_files", f) for f in os.listdir("mgf_files") if f.lower().endswith(".mgf")]
report_files = run_dia_nn(mgf_files, mzml_folder, "FirstBake_Forward", diann_path, 1.0)
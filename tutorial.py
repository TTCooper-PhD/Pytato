from pytato import *

#Step 1:
#Generate mzML files from .RAW files
raw_files="Target_Folder"
os.listdir("C:\Program Files\ProteoWizard\ProteoWizard 3.0.23")
msconvert_path="C:\Program Files\ProteoWizard\ProteoWizard 3.0.23\msconvert"
mzml_folder=convert_raw_to_mzml(raw_files,msconvert_path,output_subfolder="mzML")
#Step 2: Generate Spectral Library from 
diann_path="C:\\DIA-NN\\1.8.1\\diaNN.exe"
fasta_file="FASTA_HUMAN_Apr2023.fasta"

target_files = [os.path.join(f"{raw_files}", f) for f in os.listdir(f"{raw_files}") if f.lower().endswith(".mzml")]
library_file_list=[""]
fasta_file_list=["FASTA_HUMAN_Apr2023.fasta"]

run_dia_nn(diann_path,library_file_list,fasta_file_list,
           raw_files,"03AprTest","Test1")
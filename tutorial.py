from pytato import *

dia_nn_exe_path="C:\\DIA-NN\\1.8.1\\diaNN.exe"
msconvert_path="C:\Program Files\ProteoWizard\ProteoWizard 3.0.23\msconvert"
#Step 1:
#Generate mzML files from .RAW files
raw_file_path="Target_Folder"
msconvert_path="C:\Program Files\ProteoWizard\ProteoWizard 3.0.23\msconvert"
input_folder=convert_raw_to_mzml(raw_file_path,msconvert_path,output_subfolder="mzML")

#Step 2: First Bake(Library-free Search on Enzyme 1 followed by Targeted Search on Enzyme 2)
output_folder="Apr03_FirstBake"
report_file_name="first-bake"
fasta_file="FASTA_HUMAN_Apr2023.fasta"
library_files=[""]
fasta_files="FASTA_HUMAN_Apr2023.fasta"
cut="K*,R*"
#Enzyme1-Library-free
run_dia_nn(dia_nn_exe_path, library_files, fasta_files, input_folder, output_folder,report_file_name=report_file_name, qval=0.01,threads=30,missed_cleavages=1,
               cut="K*,R*",min_frag_mz=200,max_frag_mz=1800,min_pre_mz=300,max_pre_mz=1200, min_pep_len=7,max_pep_len=30,
               min_pre_z=1,max_pre_z=4,fasta_search=False,profiling="smart",MBR=True,fasta_speclib_annotation=False,frag_restrict_quant=True)
#Filter for "High Confidence" Protein IDs to generate a focused FASTA file for
firstbake_proteins=get_high_confidence_proteins()
firstbake_fasta=filter_fasta(fasta_file,"Test.fasta",firstbake_proteins)
#Search_Enzyme2_withEnzyme1_Library
run_dia_nn(dia_nn_exe_path, library_files, firstbake_fasta, input_folder, output_folder,report_file_name=report_file_name, qval=0.5,threads=30,missed_cleavages=2,
               cut="*F,*M,*V,*A,*I,*L",min_frag_mz=200,max_frag_mz=1800,min_pre_mz=300,max_pre_mz=1200, min_pep_len=7,max_pep_len=30,
               min_pre_z=1,max_pre_z=4,fasta_search=False,profiling="smart",MBR=True,fasta_speclib_annotation=False,frag_restrict_quant=True)


#Step 3: Second Bake

output_folder="Apr03_SecondBake"
report_file_name2="second-bake"
fasta_file="FASTA_HUMAN_Apr2023.fasta"
library_files=[""]
fasta_files="FASTA_HUMAN_Apr2023.fasta"
run_dia_nn(dia_nn_exe_path, library_files, fasta_files, input_folder, output_folder,report_file_name=report_file_name2, qval=0.01,threads=30,missed_cleavages=1,
               cut="K*,R*",min_frag_mz=200,max_frag_mz=1800,min_pre_mz=300,max_pre_mz=1200, min_pep_len=7,max_pep_len=30,ms2_acc=20,ms2_acc_cal=20,ms1_acc=20,
               min_pre_z=1,max_pre_z=4,fasta_search=False,profiling="smart",MBR=True,fasta_speclib_annotation=False,frag_restrict_quant=True,heuristic_search=True)
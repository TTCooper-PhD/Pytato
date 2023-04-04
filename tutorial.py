from pytato import *

dia_nn_exe_path="C:\\DIA-NN\\1.8.1\\diaNN.exe"
msconvert_path="C:\Program Files\ProteoWizard\ProteoWizard 3.0.23\msconvert"
#Step 1:
#Generate mzML files from .RAW files
enzyme1_file_path="Enzyme1"
msconvert_path="C:\Program Files\ProteoWizard\ProteoWizard 3.0.23\msconvert"
enzyme1_folder=convert_raw_to_mzml(enzyme1_file_path,msconvert_path,output_subfolder="mzML")
enzyme2_file_path="Enzyme2"
msconvert_path="C:\Program Files\ProteoWizard\ProteoWizard 3.0.23\msconvert"
enzyme1_folder=convert_raw_to_mzml(enzyme2_file_path,msconvert_path,output_subfolder="mzML")

#Step 2: First Bake(Library-free Search on Enzyme 1 followed by Targeted Search on Enzyme 2)
output_folder="FirstBake"
report_file_name1="trypsin"
report_file_name2="thermolysin"
fasta_file="FASTA_HUMAN_Apr2023.fasta"
library_files=[""]
fasta_files="FASTA_HUMAN_Apr2023.fasta"
cut="K*,R*"
cut2="*F,*M,*V,*A,*I,*L"

#Enzyme1_Library-Free
run_dia_nn(dia_nn_exe_path, library_files, fasta_files, enzyme1_folder, output_folder,report_file_name=report_file_name1, qval=0.01,threads=30,missed_cleavages=1,
               cut="K*,R*",min_frag_mz=200,max_frag_mz=1800,min_pre_mz=300,max_pre_mz=1200, min_pep_len=7,max_pep_len=30,ms2_acc=20,ms2_acc_cal=20,ms1_acc=20,
               min_pre_z=1,max_pre_z=4,fasta_search=True,profiling="smart",MBR=True,fasta_speclib_annotation=False,frag_restrict_quant=True,heuristic_search=True)
#Filter for "High Confidence" Protein IDs to generate a targeted FASTA file for "2nd Bake"
firstbake_proteins=get_high_confidence_proteins("Apr03_FirstBake\FirstBake_3Apr2023.tsv",fdr_threshold=0.01)
print(f"Identified {len(firstbake_proteins)} High-confidence proteins in First Bake (Enzyme1)....")
firstbake_fasta_name="Targeted_Enzyme1.fasta"
firstbake_fasta=filter_fasta(fasta_file,firstbake_fasta_name,firstbake_proteins)



#Enzyme2_Library-Free
run_dia_nn(dia_nn_exe_path, library_files, fasta_files, enzyme2_file_path, output_folder,report_file_name=report_file_name2, qval=0.01,threads=30,missed_cleavages=1,
               cut="K*,R*",min_frag_mz=200,max_frag_mz=1800,min_pre_mz=300,max_pre_mz=1200, min_pep_len=7,max_pep_len=30,ms2_acc=20,ms2_acc_cal=20,ms1_acc=20,
               min_pre_z=1,max_pre_z=4,fasta_search=True,profiling="smart",MBR=True,fasta_speclib_annotation=False,frag_restrict_quant=True,heuristic_search=True)
#Filter for "High Confidence" Protein IDs to generate a targeted FASTA file for "2nd Bake"
secondbake_proteins=get_high_confidence_proteins("Apr03_FirstBake\FirstBake_3Apr2023.tsv",fdr_threshold=0.01)
print(f"Identified {len(secondbake_proteins)} High-confidence proteins in First Bake (Enzyme2)....")
secondbake_fasta_name="Targeted_Enzyme2.fasta"
secondbake_fasta=filter_fasta(fasta_file,secondbake_fasta_name,secondbake_proteins)



#Step 3: Second Bake
report_file_name1="targeted_enzyme1"

run_dia_nn(dia_nn_exe_path, library_files, firstbake_fasta, enzyme2_file_path, output_folder,report_file_name=report_file_name1, qval=0.1,threads=30,missed_cleavages=2,
               cut="K*,R*",min_frag_mz=200,max_frag_mz=1800,min_pre_mz=300,max_pre_mz=1200, min_pep_len=7,max_pep_len=30,
               min_pre_z=1,max_pre_z=4,fasta_search=True,profiling="smart",MBR=True,fasta_speclib_annotation=False,frag_restrict_quant=True)

report_file_name2="targeted_enzyme1"
run_dia_nn(dia_nn_exe_path, library_files, firstbake_fasta, enzyme2_file_path, output_folder,report_file_name=report_file_name2, qval=0.1,threads=30,missed_cleavages=2,
               cut="K*,R*",min_frag_mz=200,max_frag_mz=1800,min_pre_mz=300,max_pre_mz=1200, min_pep_len=,max_pep_len=30,
               min_pre_z=1,max_pre_z=4,fasta_search=True,profiling="smart",MBR=True,fasta_speclib_annotation=False,frag_restrict_quant=True)
from pytato import *

#Step 1:
#Generate mzML files from .RAW files
raw_files="Target_Folder"
msconvert_path="C:\Program Files\ProteoWizard\ProteoWizard 3.0.23018.60066e9\msconvert"
mzml_folder=convert_raw_to_mzml(raw_files,msconvert_path,output_subfolder="mzML")
mzml_folder="Target_Folder\mzML"
#Step 2:
#Generate Spectral Library Using mzML files Generated in Step 1
dia_umpire_path="C:\\Users\\tcoop\\Desktop\\DIA_Umpire\\DIA_Umpire_SE.jar"
generate_spectral_library(mzml_folder,"mgf_files",dia_umpire_path)

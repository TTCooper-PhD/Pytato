
from bs4 import BeautifulSoup
import os
import requests

#``````````````````````````````````````````````````````````````````````````````````````````````````````````````


def download_msconvert(output_folder):
    url = "http://proteowizard.sourceforge.net/downloads.shtml"
    response = requests.get(url)
    response.raise_for_status()
    
    soup = BeautifulSoup(response.text, 'html.parser')
    links = soup.find_all('a')
    
    windows_download_link = None
    
    for link in links:
        href = link.get('href')
        if href and "proteowizard-windows-x86_64" in href:
            windows_download_link = href
            break
            
    if not windows_download_link:
        raise Exception("Could not find the download link for the latest ProteoWizard version for Windows.")
        
    print("Downloading ProteoWizard...")
    response = requests.get(windows_download_link)
    response.raise_for_status()
    
    zip_file = os.path.join(output_folder, "ProteoWizard.zip")
    with open(zip_file, 'wb') as f:
        f.write(response.content)
        
    print("ProteoWizard downloaded successfully. Please unzip the archive and find the msconvert.exe in the 'tools' directory.")

# # Example usage
# pytato_folder = "/path/to/Pytato"
# download_msconvert(pytato_folder)

#``````````````````````````````````````````````````````````````````````````````````````````````````````````````


def download_dia_umpire(py_path, dia_umpire_url):
    dia_umpire_jar_name = 'DIA_Umpire_SE.jar'
    dia_umpire_jar_path = os.path.join(py_path, dia_umpire_jar_name)

    if not os.path.exists(dia_umpire_jar_path):
        print(f"{dia_umpire_jar_name} not found. Downloading from GitHub...")
        response = requests.get(dia_umpire_url)
        response.raise_for_status()

        with open(dia_umpire_jar_path, 'wb') as f:
            f.write(response.content)
        print(f"{dia_umpire_jar_name} downloaded successfully.")
    else:
        print(f"{dia_umpire_jar_name} already exists in the Pytato folder.")

# Example usage
# pytato_folder_path = "/Dekstop/Pytato"
# dia_umpire_url = "https://github.com/diaumpire/DIA-Umpire/releases/download/v2.2/DIA_Umpire_SE.jar"

#````````````````````````````````````````````````````````````````````````````````````````````````````````

def download_dia_nn(output_folder, dia_nn_url):
    response = requests.get(dia_nn_url)
    response.raise_for_status()

    exe_path = os.path.join(output_folder, "dia_nn.exe")
    with open(exe_path, "wb") as f:
        f.write(response.content)
        
    print(f"DIA-NN downloaded successfully to {exe_path}")

# Example usage
# pytato_folder = "/path/to/Pytato"
# dia_nn_url = "https://github.com/vdemichev/DiaNN/releases/download/1.7.12/DiaNN-1.7.12.exe"

#``````````````````````````````````````````````````````````````````````````````````````````````````````````````


#``````````````````````````````````````````````````````````````````````````````````````````````````````````````



#``````````````````````````````````````````````````````````````````````````````````````````````````````````````


#``````````````````````````````````````````````````````````````````````````````````````````````````````````````



#``````````````````````````````````````````````````````````````````````````````````````````````````````````````



#``````````````````````````````````````````````````````````````````````````````````````````````````````````````

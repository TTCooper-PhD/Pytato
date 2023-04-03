# Pytato: A Dual-Enzyme Search Engine for DIA Proteomics

Pytato is a dual-enzyme proteomics search engine that leverages the information provided by complementary proteolysis to improve peptide and protein identification in data-independent proteomics experiments. 

## Experimental Example

Enyzme1=Trypsin
Enzyme2=Thermolysin

## Features

- Integration with DIA-NN and Percolator
- Custom Python-based pipeline for streamlined analysis
- Support for multiple proteases to improve protein identification confidence
- Generation of theoretical spectra for improved peptide matching

![logo_small](https://user-images.githubusercontent.com/36017084/229610464-03d73a08-c55e-4e9f-8dec-ac0af352a945.png)


## Installation

Before using Pytato, you need to download and set up the following external components:


1. **DIA-NN**: Download from https://github.com/vdemichev/DiaNN
2. **Percolator**: Download from https://github.com/percolator/percolator

Follow the installation instructions provided by each component's respective GitHub page or documentation.

Next, clone the Pytato repository:

```bash
git clone https://github.com/TTCooper-PhD/Pytato.git
```

Install the required python packages:

```bash
pip install -r requirements.txt
```

## Usage
- Prepare your Enzyme1 and Enzyme2-digested samples and acquire LC-MS/MS data in DIA mode.
### For a list of available enzymes and their cleavage rules, see the [Available Enzymes](enzymes.md) document.

- Process raw data files and generate spectral libraries for both Enzyme1 and Enzyme2 samples.
- Run Pytato by providing the necessary input files and parameters.

```bash
python pytato.py --enzyme1_data trypsin_data.mzML --enzyme2_data thermolysin_data.mzML  --output output_directory
```

For detailed usage instructions and available options, refer to the documentation.

License
Pytato is released under the Apache License, Version 2.0. Please note that the external components (MS-Fragger, DIA-Umpire, DIA-NN, and Percolator) have their own licenses, which you should review before using them in conjunction with Pytato.

Acknowledgements
We would like to acknowledge the developers of MS-Fragger, DIA-Umpire, DIA-NN, and Percolator for their valuable contributions to the field of proteomics. Please refer to their respective GitHub repositories and publications for more information.

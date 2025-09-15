# Repository Structure


Interaction_data_updated.zip   # contains input data for conformers script to analyse interaction energy(N1, N2, P1, P2)  
Internal_coordinates.py        # main script to extract and analys the internal coordinates(bond angle,bond length, dihedral angle)
README.md                      # project overview and background  
Tool Usage.md                  # repository structure and usage instructions 
  

## Usage

    Clone the repository:

git clone https://github.com/your-username/MASTER-S-PROJECT.git
cd MASTER-S-PROJECT

Unzip the data:

unzip Interaction_data_updated.zip

Run the analysis script:

    python Interaction_ver_2.py

Input files required (per conformer folder N1, N2, P1, P2):

    .log file (quantum chemistry output)

    .pdb file (coordinates)

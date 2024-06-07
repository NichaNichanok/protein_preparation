import os
import numpy as np
import pymol
import __main__

def get_pdb(pdb_code=""):
    """
    Download or load the protein structure PDB file.

    **Args:**
        pdb_code (string): PDB code according to RCSB (e.g., 6o0k) or the AlphaFold database (e.g., AF-Q16611-F1-model_v4), 
                           or leave as an empty string and provide the path to the PDB file.

    **Returns:**
        pdb_file (string): Path to the protein structure PDB file.
    """
    if pdb_code is None or pdb_code == "":
        pdb_file = input("Please provide the path to your PDB file: ")
        return pdb_file
    elif os.path.isfile(pdb_code):
        return pdb_code
    elif len(pdb_code) == 4:
        os.system(f"wget -qnc https://files.rcsb.org/view/{pdb_code}.pdb")
        print(f"Downloaded the {pdb_code} PDB structure from RCSB website...")
        return f"{pdb_code}.pdb"
    else:
        os.system(f"wget -qnc https://alphafold.ebi.ac.uk/files/AF-{pdb_code}-F1-model_v4.pdb")
        print(f"Downloaded the {pdb_code} PDB structure from AlphaFold database website...")
        return f"AF-{pdb_code}-F1-model_v4.pdb"

def grid_size(target_pdb, pymol_cmd, size=34):
    """
    Calculate the grid coordinate (center of mass) from the list of selected residues.

    **Args:**
        target_pdb (string): PDB code (according to RCSB or the AlphaFold database), 
                             or the path to the PDB file.
        pymol_cmd (string): Selected list of binding residues, e.g., "resi 112 + resi 137 + resi 149 + resi 115".

    **Returns:**
        size(int): grid size, calculate by diameter of the selected binding residues + the distance between the centerofmass of protein and the selected binding residues

    """
    protein_structure = get_pdb(target_pdb)
    
    __main__.pymol_argv = ['pymol', '-qc']  # Quiet and no GUI
    pymol.finish_launching()

    # Load protein structure
    pymol.cmd.load(protein_structure, 'target_protein')

    # Compute the overall center of mass
    binding_res_coords = pymol.cmd.centerofmass(pymol_cmd)
    protein_coords = pymol.cmd.centerofmass('target_protein')
    
    # Print the overall center of mass
    print("Overall center of mass of binding residues:", binding_res_coords)
    print("Overall center of mass of protein:", protein_coords)

    # Calculate the Euclidean distance between the binding residues and the protein center of mass
    bres_point = np.array(binding_res_coords)
    prot_point = np.array(protein_coords)
    euc_distance = np.linalg.norm(bres_point - prot_point)
    print(f"Distance between the binding residues and protein = {euc_distance} nm")


    #get diameter of protein
    # Get coordinates of all atoms in the protein
    atom_coords = pymol.cmd.get_coords("target_protein")
    distances = np.linalg.norm(atom_coords[:, None] - atom_coords, axis=-1)
    np.fill_diagonal(distances, 0)    # Set the diagonal elements to 0 to avoid self-distances
    diameter = np.max(distances)
    print("The diameter of the protein is:", diameter)

    #get diameter of binding residues
    pymol.cmd.select('binding_res', pymol_cmd)
    # Get coordinates of all atoms in the protein
    atom_coords_bres = pymol.cmd.get_coords("binding_res")
    # Calculate the pairwise distances between all atoms
    distances_bres = np.linalg.norm(atom_coords_bres[:, None] - atom_coords_bres, axis=-1)
    np.fill_diagonal(distances_bres, 0)# Set the diagonal elements to 0 to avoid self-distances
    diameter_bres = np.max(distances_bres)
    # Print the diameter
    print("The diameter of the binding residues is:", diameter_bres)

    size = round(diameter_bres + euc_distance)

    # Quit PyMOL
    pymol.cmd.quit()

    return size


if __name__ == "__main__":
    pdb_file = get_pdb("6o0k")
    binding_res_selection = "resi 112 + resi 137 + resi 149 + resi 115 + resi 133 + resi 108 + resi 104 + resi 152 + resi 153 + resi 146 + resi 111 + resi 156 + resi 136 + resi 145 + resi 148"
    
    grid_size = grid_size("6o0k", binding_res_selection)
    print(grid_size)


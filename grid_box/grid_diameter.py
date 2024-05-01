import pymol

def grid_diameter(input_path, ligand_name):
    """
    Calculate the diameter: maximum distance between any two atoms in the specified ligand.

    This function initializes PyMOL, loads the specified structure file, selects the specified ligand excluding protein atoms,
    retrieves the list of atom names in the ligand, calculates the maximum distance between any two atoms, and prints
    the result.

    Parameters:
    - input_path (str): The file path of the structure file.
    - ligand_name (str): The name of the ligand in the structure file.

    Output: diameter in Angstrom, which will be used for the grid size
    """
    # Initialize PyMOL
    pymol.finish_launching(['pymol', '-qc'])

    # Load the specified structure file
    pymol.cmd.load(input_path)

    # Select the specified ligand excluding protein atoms
    pymol.cmd.select(ligand_name, 'hetatm')

    # Get the list of atom names in the specified ligand
    atom_names = [atom.name for atom in pymol.cmd.get_model(ligand_name).atom]

    # Initialize variable to store maximum distance
    max_distance = 0

    # Iterate through each pair of atoms
    for i, atom_name1 in enumerate(atom_names):
        for atom_name2 in atom_names[i+1:]:
            # Calculate distance between the two atoms
            distance = pymol.cmd.distance(f'{atom_name1}_to_{atom_name2}', f'resname LBM and name {atom_name1}', f'resname LBM and name {atom_name2}')
            # Update maximum distance if needed
            print(f"{i} distance: {distance}")
            max_distance = max(max_distance, distance)
            print(f"max distance: {max_distance}")
    # Quit PyMOL
    pymol.cmd.quit()

    return round(max_distance)


    

if __name__ == "__main__":
    input_path = '/home/nauevech/Documents/protein_preparation/protein_preparation/protein_preprocessing/input/6o0k.pdb'
    ligand_name = 'dummie'  
    diameter = grid_diameter(input_path, ligand_name)
    print("Diameter of the ligand:", diameter, "Angstroms")
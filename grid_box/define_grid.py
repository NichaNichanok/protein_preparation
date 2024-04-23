def define_grid(input_path, output_directory, bybinding_res = True, byligand = False):
    """
    Process experimental holo-PDB files in the input directory:
    - identify the ligand and it center of mass for further use as grid coordinate (config.txt)
    - save config file with the coordinate in the output directory.
    
    Args:
    - input_path (str): Path to the directory containing PDB files.
    - output_directory (str): Path to the directory for saving modified PDB files.
    """
    import __main__
    __main__.pymol_argv = ['pymol', '-qc']  # Quiet and no GUI
    # Initialize PyMOL
    pymol.finish_launching()

    # Loop over all files in the input directory
    for filename in os.listdir(input_path):
        if filename.endswith(".pdb"):
            try:
                # Load the PDB file
                pdb_file_path = os.path.join(input_path, filename)
                # Remove the .pdb extension to create a valid object name
                object_name = os.path.splitext(filename)[0]
                pymol.cmd.load(pdb_file_path, object_name)

                # Select the ligand and calculate its center of mass
                pymol.cmd.select("ligand", "organic")
                center_of_mass = pymol.cmd.centerofmass("ligand")

                # Save the coordinates of the ligand's center of mass to a text file
                output_config_path = os.path.join(output_directory, "config.txt")
                with open(output_config_path, "w") as config_file:
                    config_file.write(f"The grid coordinates of '{filename.split('.')[0]}' protein by its true ligand:\n")
                    config_file.write("X: {:.3f}\n".format(center_of_mass[0]))
                    config_file.write("Y: {:.3f}\n".format(center_of_mass[1]))
                    config_file.write("Z: {:.3f}\n".format(center_of_mass[2]))

                print(f"The grid coordinates of '{filename}' protein by its true ligand:", center_of_mass)

                print(f"Processed {filename}. Output saved to {output_config_path}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")

    # Quit PyMOL
    pymol.cmd.quit()

def define_grid_bybindingres(input_path, csv_path, output_directory, pbind=0.8):
    """
    Define the grid from the predicted binding site:
    - Select the binding residues (as binding_res) based on the pbind value by default 0.8
    - Calculate the center of mass of the binding_res for further use as grid coordinate (config.txt)
    - Save config file with the coordinate in the output directory.
    
    Args:
    - input_path (str): Path to the directory containing PDB files.
    - csv_path (str): Path to the predicted binding site by af2bind in CSV format.
    - output_directory (str): Path to the directory for saving modified PDB files.
    - pbind (foat): pbind value cutoff
    """
    __main__.pymol_argv = ['pymol', '-qc']  # Quiet and no GUI
    pymol.finish_launching()

    # Load protein structure
    pymol.cmd.load(input_path, 'protein')

    # Initialize list to store binding residues
    binding_res = set()

    # Read CSV file and extract residues with p(bind) > 0.8
    with open(csv_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            data = line.strip().split(',')
            resi = int(data[2])  # Convert to integer
            p_bind = float(data[4])  
            if p_bind > pbind:
                binding_res.add(resi)  # Add residue to binding_res set

    # Select all binding residues
    pymol.cmd.select('binding_res', 'resi ' + '+'.join(map(str, binding_res)))
    #pymol.cmd.show('sphere', 'binding_res')
    #pymol.cmd.show('sticks', 'binding_res')  # Show sticks for binding residues
    #pymol.cmd.color('red', 'binding_res')

    # Compute the overall center of mass
    binding_res_coords = pymol.cmd.centerofmass('resi ' + '+'.join(map(str, binding_res)))


    # Print the overall center of mass
    print(f"The grid coordinates of '{protein_name}' protein by selected binding residues with the pbind > {pbind}:", binding_res_coords)

    protein_name= input_path.split("/")[-1].split(".")[0]

    # Save the coordinates of the binding residues' center of mass to a text file
    output_config_path = os.path.join(output_directory, "config.txt")
    with open(output_config_path, "w") as config_file:
        config_file.write(f"The grid coordinates of '{protein_name}' protein by selected binding residues with the pbind > {pbind}:\n")
        config_file.write("X: {:.3f}\n".format(binding_res_coords[0]))
        config_file.write("Y: {:.3f}\n".format(binding_res_coords[1]))
        config_file.write("Z: {:.3f}\n".format(binding_res_coords[2]))

    print(f"Output saved to {output_config_path}")

    # Quit PyMOL
    pymol.cmd.quit()
import pymol 
import os
import __main__


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
    pymol.cmd.show('sphere', 'binding_res')
    pymol.cmd.show('sticks', 'binding_res')  # Show sticks for binding residues
    pymol.cmd.color('red', 'binding_res')

    # Compute the overall center of mass
    binding_res_coords = pymol.cmd.centerofmass('resi ' + '+'.join(map(str, binding_res)))


    # Print the overall center of mass
    print("Overall center of mass of binding residues:", binding_res_coords)

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


def visualize_binding_residues(input_path, csv_path, pbind=0.8):
    pymol.finish_launching()

    # Load protein structure
    pymol.cmd.load(input_path, 'protein')
    pymol.cmd.hide()
    pymol.cmd.bg_color("white")
    pymol.cmd.show("mesh")
    pymol.cmd.color("blue")

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
    pymol.cmd.show('sphere', 'binding_res')
    pymol.cmd.show('sticks', 'binding_res')  # Show sticks for binding residues
    pymol.cmd.color('red', 'binding_res')

    # Compute the overall center of mass
    binding_res_coords = pymol.cmd.centerofmass('resi ' + '+'.join(map(str, binding_res)))

    # Print the overall center of mass
    print("Overall center of mass of binding residues:", binding_res_coords)

if __name__ == "__main__":
    input_path = '/home/nauevech/Documents/protein_preparation/6o0k.pdb'
    csv_path = '/home/nauevech/Documents/protein_preparation/af2bind_6o0k_results.csv'
    output_directory = '/home/nauevech/Documents/protein_preparation/'
    #visualize_binding_residues(input_path, csv_path)
    define_grid_bybindingres(input_path, csv_path, output_directory)

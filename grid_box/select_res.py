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


def visualize_binding_residues(input_path, csv_path, pbind=0.8):
    pymol.finish_launching()

    # Load protein structure
    pymol.cmd.load(input_path, 'protein')
    pymol.cmd.hide()
    pymol.cmd.bg_color("black")
    pymol.cmd.show("surface")
    pymol.cmd.color("white")
    pymol.cmd.select('LBM', 'resn LBM')
    pymol.cmd.show('sticks', 'LBM')
    pymol.cmd.util.cba('orange', 'LBM')


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
    pymol.cmd.select(f'binding_res{pbind}', 'resi ' + '+'.join(map(str, binding_res)))
    pymol.cmd.show('sphere', f'binding_res{pbind}')
    pymol.cmd.show('sticks', f'binding_res{pbind}')  # Show sticks for binding residues
    pymol.cmd.color('red', f'binding_res{pbind}')

    # Compute the overall center of mass
    binding_res_coords = pymol.cmd.centerofmass('resi ' + '+'.join(map(str, binding_res)))

    # Print the overall center of mass
    print("Overall center of mass of binding residues:", binding_res_coords)

    # Define grid coordinates based on center of mass and size
    center_x, center_y, center_z = binding_res_coords
    size = 34  # Size of the grid box

    # Calculate vertices of the box
    x_min = center_x - (size / 2)
    x_max = center_x + (size / 2)
    y_min = center_y - (size / 2)
    y_max = center_y + (size / 2)
    z_min = center_z - (size / 2)
    z_max = center_z + (size / 2)

    # Create a box object
    box_name = "grid_box"
    pymol.cmd.pseudoatom(box_name, pos=[0, 0, 0])  # Create a dummy atom to serve as the center
    pymol.cmd.create(box_name, "byres ((x >= {}) and (x <= {})) and ((y >= {}) and (y <= {})) and ((z >= {}) and (z <= {}))".format(x_min, x_max, y_min, y_max, z_min, z_max))

    # Show the box
    pymol.cmd.show_as("surface", box_name)

    # Center and zoom the view
    pymol.cmd.zoom(box_name)

    # Color the box
    pymol.cmd.color("gray", box_name)

    # Ensure the box is visible
    pymol.cmd.show()

if __name__ == "__main__":
    input_path = '/home/nauevech/Documents/protein_preparation/6o0k.pdb'
    csv_path = '/home/nauevech/Documents/protein_preparation/af2bind_6o0k_results.csv'
    output_directory = '/home/nauevech/Documents/protein_preparation/'
    visualize_binding_residues(input_path, csv_path, 0.8)
    #define_grid_bybindingres(input_path, csv_path, output_directory)

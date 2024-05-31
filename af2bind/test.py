import pymol 
import os
import __main__



def grid_coordinate(target_pdb, selected_resi, size=34):
    
    __main__.pymol_argv = ['pymol', '-qc']  # Quiet and no GUI
    pymol.finish_launching()

    # Load protein structure
    pymol.cmd.load(target_pdb, 'protein')

    # Initialize list to store binding residues
    binding_res = set()
    # Compute the overall center of mass
    binding_res_coords = pymol.cmd.centerofmass(selected_resi)

    # Print the overall center of mass
    print("Overall center of mass of binding residues:", binding_res_coords)

    # Define grid coordinates based on center of mass and size
    center_x, center_y, center_z = binding_res_coords
    size = 34  # Size of the grid box
    
    print(f"center_x: {center_x}")
    print(f"center_y: {center_y}")
    print(f"center_z: {center_z}")
    
    # Print the overall center of mass
    print(f"The grid coordinates of '{target_pdb}' protein by top 15 amino acids:", binding_res_coords)

    protein_name= target_pdb.split("/")[-1].split(".")[0]

    # Save the coordinates of the binding residues' center of mass to a text file
    output_config_path = os.path.join(os.getcwd(), "config.txt")
    with open(output_config_path, "w") as config_file:
        config_file.write(f"The grid coordinates of '{target_pdb}' protein by top 15 amino acids:\n")
        config_file.write("center_x = {:.2f}\n".format(binding_res_coords[0]))
        config_file.write("center_y: {:.2f}\n".format(binding_res_coords[1]))
        config_file.write("center_z: {:.2f}\n".format(binding_res_coords[2]))

    print(f"Output saved to {output_config_path}")

    # Quit PyMOL
    pymol.cmd.quit()

#test
selected_resi = "resi 112 + resi 137 + resi 149 + resi 115 + resi 133 + resi 153 + resi 104 + resi 108 + resi 152 + resi 146 + resi 111 + resi 156 + resi 136 + resi 145 + resi 148"
target_pdb = "/home/nauevech/Documents/protein_preparation/protein_preparation/af2bind/6o0k.pdb"
grid_coordinate(target_pdb,selected_resi)
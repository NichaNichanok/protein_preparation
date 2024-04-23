import os
import pymol
import __main__

# This function only works if there is only ONE ligand in the using holo protein structure!

def define_grid_byligand(input_path, output_directory):
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




def define_grid_rmnp(input_path, output_directory):
    """
    Process experimental PDB files in the input directory:
    - identify the ligand and it center of mass for further use as grid coordinate (config.txt)
    - remove non-protein atoms
    - save modified files in the output directory.
    
    Args:
    - input_path (str): Path to the directory containing PDB files.
    - output_directory (str): Path to the directory for saving modified PDB files.
    """

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

                # Remove non-protein atoms
                pymol.cmd.remove("not polymer.protein")

                # Save the modified structure
                output_filename = os.path.splitext(filename)[0] + "_rmnpn.pdb"
                output_file_path = os.path.join(output_directory, output_filename)
                pymol.cmd.save(output_file_path, object_name, format="pdb")

                print(f"Processed {filename}. Output saved to {output_file_path}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")

    # Quit PyMOL
    pymol.cmd.quit()


if __name__ == "__main__":
    input_path = "/home/nauevech/Documents/protein_preparation/input_pdb_files/oneligand"
    output_directory = "/home/nauevech/Documents/protein_preparation/output_pdb_files/"
    define_grid_byligand(input_path, output_directory)
import os, pymol
import __main__

def pdb_processed(input_path, output_directory, pH = 7.4):
    """
    Process experimental PDB files in the input directory:
    - Identify the ligand and its center of mass for further use as grid coordinate
    - Remove non-protein atoms, protonate at pH 7.4 and convert to pdbqt
    - Save the ligand file in PDBQT format
    - The grid coordinate and the smile format are save in config.txt
    - Save modified files in the output directory.
    
    Args:
    - input_path (str): Path to the directory containing PDB files.
    - output_directory (str): Path to the directory for saving modified PDB files.
    - protonate pH, 7.4 by default.
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
                    config_file.write(f"Protein PDB-ID: {filename[:-4]}\n")
                    config_file.write("Grid box coordinates by center of mass of its true ligand:\n")
                    config_file.write("X: {:.3f}\n".format(center_of_mass[0]))
                    config_file.write("Y: {:.3f}\n".format(center_of_mass[1]))
                    config_file.write("Z: {:.3f}\n".format(center_of_mass[2]))

                    # Convert the ligand to PDBQT format
                    ligand_pdb = os.path.join(output_directory, "ligand.pdb")
                    pymol.cmd.save(ligand_pdb, "ligand", format="pdb")
                    output_pdbqt = os.path.join(output_directory, "ligand.pdbqt")
                    os.system(f"obabel {ligand_pdb} -opdbqt -xr -O {output_pdbqt}")

                    # Obtain the SMILES representation of the ligand
                    output_smi = os.path.join(output_directory, "ligand.smi")
                    os.system(f"obabel {ligand_pdb} -osmi -O {output_smi} --gen3d -h")
                    with open(output_smi, "r") as smi_file:
                        smiles = smi_file.read().strip()
                    config_file.write("SMILES: {}\n".format(smiles))

                # Remove non-protein atoms
                pymol.cmd.remove("not polymer.protein")

                # Save the pdb structure without non-protein molecules
                output_filename_rmnpm = os.path.splitext(filename)[0] + "_rmnpm.pdb"
                output_file_path_rmnpm = os.path.join(output_directory, output_filename_rmnpm)
                pymol.cmd.save(output_file_path_rmnpm, format="pdb")

                output_filename_protonated = os.path.splitext(filename)[0] + "_protonated.pdb"
                output_file_path_protonated = os.path.join(output_directory, output_filename_protonated)
                # Run the command to protonate at pH 7.4, save in pdbqt
                os.system(f"obabel {output_file_path_rmnpm} -opdb -h -p {pH} -O {output_file_path_protonated}")

                # Save the final processed file in pdbqt format 
                output_filename_processed = os.path.splitext(filename)[0] + "_processed.pdbqt"
                output_file_path_processed = os.path.join(output_directory, output_filename_processed)
                os.system(f"obabel {output_file_path_protonated} -opdbqt -xr -O {output_file_path_processed}")

                print(f"Processed {filename}. Output saved to {output_file_path_processed}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")

    # Quit PyMOL
    pymol.cmd.quit()


if __name__ == "__main__":
    input_path = "/home/nauevech/Documents/protein_preparation/input_pdb_files/oneligand"
    output_directory = "/home/nauevech/Documents/protein_preparation/output_pdb_files/"
    pdb_processed(input_path, output_directory)

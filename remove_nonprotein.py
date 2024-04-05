import os
import pymol

def remove_nonprotein(input_path, output_directory):
    """
    Process experimental PDB files in the input directory, 
    remove non-protein atoms, and save modified files in the output directory.
    
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
    input_path = "/home/nauevech/Documents/protein_preparation/input_pdb_files"
    output_directory = "/home/nauevech/Documents/protein_preparation/output_pdb_files/"
    remove_nonprotein(input_path, output_directory)

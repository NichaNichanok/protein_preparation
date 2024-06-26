import os
# this code require the Openbabel install

def energy_minimize_pdbqt(input_directory, output_directory):
    """
    Process PDB files in the input directory, 
    perform the energy minimize, and convert them to PDBQT format.
    
    Args:
    - input_directory (str): Path to the directory containing PDB files.
    - output_directory (str): Path to the directory for saving modified PDBQT files.
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Loop over all files in the input directory
    for filename in os.listdir(input_directory):
        # Check if the file is a PDB file
        if filename.endswith(".pdb"):
            try:
                # Input PDB file path
                input_pdb = os.path.join(input_directory, filename)

                # Minimized PDB file path
                minimized_pdb = os.path.join(output_directory, f"{os.path.splitext(filename)[0]}_minimized.pdb")

                # Output PDBQT file path
                output_pdbqt = os.path.join(output_directory, f"{os.path.splitext(filename)[0]}_protein.pdbqt")

                # Run the command to protonate at pH 7.4 and conduct local energy minimizer the input PDB file
                os.system(f"obabel {input_pdb} -opdb --minimize  -O {minimized_pdb}")

                # Run the command to convert PDB file to PDBQT format
                os.system(f"obabel {minimized_pdb} -opdbqt -xr -O {output_pdbqt}")
                
                print(f"Processed {filename}. minimized file saved to {minimized_pdb}. PDBQT file saved to {output_pdbqt}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")

if __name__ == "__main__":
    # Get the current working directory
    current_directory = os.getcwd()

    # Input directory containing PDB files
    input_directory = os.path.join(current_directory, "output_pdb_files")

    # Output directory for saving modified files
    output_directory = os.path.join(current_directory, "output_files")

    energy_minimize_pdbqt(input_directory, output_directory)

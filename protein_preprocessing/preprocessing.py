import os, pymol
import __main__
import requests
from bs4 import BeautifulSoup
#from fetch_rcsb import fetch_ligand_name

def fetch_ligand_name(pdb_id):
    """
    Search for the drug-ligand name defined in the PDB file 
    by fetching the information in the website, 
    i.e. the binding affinity annotations ID
    
    Args:
    - PDB ID: 4 letter code protein from the RCSB.org website
    
    Output:
    ligand_name (str) which defined in the pdb format
    return none: if there is no further binding affinity annotation, i.e. apo form
    """
    # URL for the RCSB structure page with the provided PDB ID
    url = f"https://www.rcsb.org/structure/{pdb_id}"

    # Send a GET request to the URL
    response = requests.get(url)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Parse the HTML content of the page
        soup = BeautifulSoup(response.content, "html.parser")

        # Find the table containing the binding affinity annotations
        table = soup.find("table", {"class": "table table-bordered table-condensed", "id": "binding-affinity-table"})

        # Check if the table was found
        if table:
            # Find the first row in the table body
            row = table.find("tbody").find("tr")

            # Find the cell containing the ID
            cell = row.find("td")

            # Extract the text from the cell
            id_name = cell.get_text(strip=True)

            return id_name
        else:
            return None

def crystal_processing(input_path, output_directory, pH = 7.4):
    """
    Process experimental PDB files (in holo format) in the input directory:
    - Identify the ligand and its center of mass for further use as grid coordinate
    - Remove non-protein atoms, protonate at pH 7.4 and convert to pdbqt
    - Save the ligand file in PDBQT format
    - The grid coordinate and the smile format are save in config.txt
    - Save modified files in the output directory.
    
    Args:
    - input_path (str): Path to the directory containing PDB files.
                        * filename have to be a PDB-id e.g. 6o0k_example
    - output_directory (str): Path to the directory for saving modified PDB files.
    - ligand_name (str): uppercase of 3-letter ligand name, used in RCSB
    - protonate pH, 7.4 by default.
    """
    __main__.pymol_argv = ['pymol', '-qc']
    pymol.finish_launching()

    # Loop over all files in the input directory
    for filename in os.listdir(input_path):
        if filename.endswith(".pdb"):
            try:
                # Load the PDB file
                pdb_file_path = os.path.join(input_path, filename)
                object_name = os.path.splitext(filename)[0]
                pymol.cmd.load(pdb_file_path, object_name)

                # Fetch the ligand name
                ligand = fetch_ligand_name(filename[0:4])
                pymol.cmd.select("ligand", f"resn {ligand}")
                center_of_mass = pymol.cmd.centerofmass("ligand")

                # Calculate the ligand diameter and grid size
                atom_names = [atom.name for atom in pymol.cmd.get_model("ligand").atom]
                max_distance = 0
                for i, atom_name1 in enumerate(atom_names):
                    for atom_name2 in atom_names[i+1:]:
                        distance = pymol.cmd.distance(f'{atom_name1}_to_{atom_name2}',
                                                      f'resname {ligand} and name {atom_name1}',
                                                      f'resname {ligand} and name {atom_name2}')
                        max_distance = max(max_distance, distance)
                diameter = round(max_distance)
                size = round(16 + 0.8*diameter)

                ## Process ligand
                ligand_filename_pdb = os.path.splitext(filename)[0] + "_ligand.pdb"
                ligand_pdb = os.path.join(output_directory, ligand_filename_pdb)
                pymol.cmd.save(ligand_pdb, "ligand", format="pdb")
                ligand_filename_pdbqt = os.path.splitext(filename)[0] + "_ligand.pdbqt"
                output_pdbqt = os.path.join(output_directory, ligand_filename_pdbqt)
                ligand_filename_smi = os.path.splitext(filename)[0] + "_ligand.smi"
                output_smi = os.path.join(output_directory, ligand_filename_smi)
                # add gastaiger charges, set TORDOF, convert to pdbqt and smi
                os.system(f"obabel {ligand_pdb} -O {output_pdbqt} --gen3d -p TORDOF --partialcharge gasteiger")

                ## Process protein
                # remove non-protein molecules
                pymol.cmd.remove("not polymer.protein")
                output_filename_rmnpm = os.path.splitext(filename)[0] + "_rmnpm.pdb"
                output_file_path_rmnpm = os.path.join(output_directory, output_filename_rmnpm)
                pymol.cmd.save(output_file_path_rmnpm, format="pdb")
                # protonate at pH 7.4 by default,use for partial charges (eem is Bultnck B3LYP/6-13G*/MPA)
                output_filename_protonated = os.path.splitext(filename)[0] + "_protonated.pdb"
                output_file_path_protonated = os.path.join(output_directory, output_filename_protonated)
                os.system(f"obabel {output_file_path_rmnpm} -opdb -p {pH} -partialcharge eem -O {output_file_path_protonated}")
                # convert to pdbqt file
                output_filename_processed = os.path.splitext(filename)[0] + "_processed.pdbqt"
                output_file_path_processed = os.path.join(output_directory, output_filename_processed)
                os.system(f"obabel {output_file_path_protonated} -opdbqt -xr -O {output_file_path_processed}")

                # Write config.txt file
                config_filename = os.path.splitext(filename)[0] + "_config.txt"
                output_config_path = os.path.join(output_directory, config_filename)
                with open(output_config_path, "w") as config_file:
                    config_file.write(f"receptor = {filename[:-4]}.pdbqt\n")
                    config_file.write(f"ligand = ligand.pdbqt\n")
                    config_file.write("center_x = {:.3f}\n".format(center_of_mass[0]))
                    config_file.write("center_y = {:.3f}\n".format(center_of_mass[1]))
                    config_file.write("center_z = {:.3f}\n".format(center_of_mass[2]))
                    config_file.write("\n")
                    config_file.write("size_x = {:.3f}\n".format(size))
                    config_file.write("size_y = {:.3f}\n".format(size))
                    config_file.write("size_z = {:.3f}\n".format(size))

                print(f"Processed {filename}. Output saved to {output_file_path_processed}")
            except Exception as e:
                print(f"Error processing {filename}: {e}")

    # Quit PyMOL
    pymol.cmd.quit()

if __name__ == "__main__":
    input_path = "/home/nauevech/Documents/protein_preparation/protein_preparation/protein_preprocessing/input/"
    output_directory = "/home/nauevech/Documents/protein_preparation/protein_preparation/protein_preprocessing/output/"
    crystal_processing(input_path, output_directory)
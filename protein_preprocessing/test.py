import streamlit as st
import pymol
import py3Dmol
import os
import requests
from bs4 import BeautifulSoup

st.sidebar.title('Protein Preparation for Virtual Screening')
st.sidebar.write('The code of this page is shown in [GitHub](https://github.com/NichaNichanok) in "pdb_fetch_app.py".')

def fetch_ligand_name(pdb_id):
    """
    Search for the ligand name defined in the PDB file 
    by fetching the information from the RCSB.org website.
    
    Args:
    - pdb_id (str): 4 letter code protein from the RCSB.org website
    
    Returns:
    - ligand_name (str): Ligand name defined in the pdb format
    - None: if there is no further binding affinity annotation, i.e. apo form
    """
    url = f"https://www.rcsb.org/structure/{pdb_id}"
    response = requests.get(url)

    if response.status_code == 200:
        soup = BeautifulSoup(response.content, "html.parser")
        table = soup.find("table", {"class": "table table-bordered table-condensed", "id": "binding-affinity-table"})

        if table:
            row = table.find("tbody").find("tr")
            cell = row.find("td")
            ligand_name = cell.get_text(strip=True)
            return ligand_name

    return None

def render_mol(pdb):
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb, 'pdb')
    
    pdbview.setStyle({'chain': 'A'}, {'cartoon': {'color': 'green'}})
    pdbview.setStyle({'resn': ['HOH', 'WAT']}, {'cross': {}})
    pdbview.setStyle({'hetflag': True}, {'stick': {'colorscheme':'yellowCarbon'}})
    pdbview.addStyle({'within': {'distance': '5', 'sel': {'organic': True}}}, {'stick': {}})
    pdbview.addSurface(py3Dmol.VDW, {'opacity': 0.85, 'color': 'white'}, {'not': {'hetflag': True}})
    
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.zoom(0.8)  
    pdbview.spin(True, speed=0.025)
    
    st.subheader('Visualization of the Protein Structure')
    st.write(pdbview)

DEFAULT_SEQ = "1sqt"
txt = st.sidebar.text_area('Input PDB-ID: with ONLY a ligand', DEFAULT_SEQ, height=275)

def get_pdb(pdb_id):
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(pdb_url)

    if response.status_code == 200:
        pdb_filename = f"{pdb_id}.pdb"
        with open(pdb_filename, 'wb') as f:
            f.write(response.content)
        return pdb_filename
    else:
        raise ValueError(f"Failed to download PDB file for {pdb_id}")

def crystal_processing(pdb_id, output_directory=os.getcwd(), pH=7.4):
    pdb_filename = get_pdb(pdb_id)
    pymol.cmd.load(pdb_filename)
    
    ligand_name = fetch_ligand_name(pdb_id)
    pymol.cmd.select("ligand", f"resn {ligand_name}")
    center_of_mass = pymol.cmd.centerofmass("ligand")

    atom_names = [atom.name for atom in pymol.cmd.get_model("ligand").atom]
    max_distance = 0
    
    for i, atom_name1 in enumerate(atom_names):
        for atom_name2 in atom_names[i + 1:]:
            distance = pymol.cmd.distance(f'{atom_name1}_to_{atom_name2}', 
                                          f'resname {ligand_name} and name {atom_name1}', 
                                          f'resname {ligand_name} and name {atom_name2}')
            max_distance = max(max_distance, distance)

    diameter = round(max_distance)
    size = round(16 + 0.8 * diameter)

    ligand_pdb = os.path.join(output_directory, "ligand.pdb")
    pymol.cmd.save(ligand_pdb, "ligand", format="pdb")
    
    output_pdbqt = os.path.join(output_directory, "ligand.pdbqt")
    os.system(f"obabel {ligand_pdb} -opdbqt -xr -O {output_pdbqt}")
    
    output_filename_rmnpm = os.path.splitext(pdb_filename)[0] + "_rmnpm.pdb"
    output_file_path_rmnpm = os.path.join(output_directory, output_filename_rmnpm)
    pymol.cmd.save(output_file_path_rmnpm, format="pdb")
    
    output_filename_protonated = os.path.splitext(pdb_filename)[0] + "_protonated.pdb"
    output_file_path_protonated = os.path.join(output_directory, output_filename_protonated)
    os.system(f"obabel {output_file_path_rmnpm} -opdb -p {pH} -h -O {output_file_path_protonated}")
    
    output_filename_processed = os.path.splitext(pdb_filename)[0] + "_processed.pdbqt"
    output_file_path_processed = os.path.join(output_directory, output_filename_processed)
    os.system(f"obabel {output_file_path_protonated} -opdbqt -xr -O {output_file_path_processed}")
    
    output_config_path = os.path.join(output_directory, "config.txt")
    with open(output_config_path, "w") as config_file:
        config_file.write(f"receptor = {pdb_filename[:-4]}.pdbqt\n")
        config_file.write(f"ligand = ligand.pdbqt\n")
        config_file.write("center_x = {:.3f}\n".format(center_of_mass[0]))
        config_file.write("center_y = {:.3f}\n".format(center_of_mass[1]))
        config_file.write("center_z = {:.3f}\n".format(center_of_mass[2]))
        config_file.write("\n")
        config_file.write("size_x = {:.3f}\n".format(size))
        config_file.write("size_y = {:.3f}\n".format(size))
        config_file.write("size_z = {:.3f}\n".format(size))

    print(f"Processed {pdb_filename}. Output saved to {output_file_path_processed}")
    pymol.cmd.quit()

def update(sequence=txt):
    pdb_id = sequence[:4].upper()
    pdb_filename = get_pdb(pdb_id)

    with open(pdb_filename, 'r') as f:
        pdb_string = f.read()

    st.subheader(f'Visualization of the {pdb_id} protein structure')
    render_mol(pdb_string)

predict = st.sidebar.button("Let's have a look!", on_click=update)
preparation = st.sidebar.button("Let's get the protein prepare", on_click=crystal_processing(pdb_id=txt[:4]))

if not predict:
    st.warning('👈 Enter protein PDB-ID!')

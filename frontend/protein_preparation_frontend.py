# This is app is inspired by Chanin Nantasenamat (Data Professor) https://youtube.com/dataprofessor


import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio

#st.set_page_config(layout = 'wide')
st.sidebar.title('Protein Preparation for virtual scrreening')
st.sidebar.write('The code of this page is shown in [*Github*](https://github.com/NichaNichanok) in "pdb_fetch_app.py.')


# stmol
def render_mol(pdb):
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb, 'pdb')
    
    # Set style for protein atoms (displayed as surface with gray color)
    pdbview.setStyle({'chain': 'A'}, {'cartoon': {'color': 'green'}})

    # Set style for water molecules (displayed as a cross)
    pdbview.setStyle({'resn': ['HOH', 'WAT']}, {'cross': {}})
    
    # Set style for HETATM atoms (displayed in stick style)
    pdbview.setStyle({'hetflag': True}, {'stick': {'colorscheme':'yellowCarbon'}})

    # Set style for atoms within 5 angstrom of ligand UH7 (displayed in stick style)
    pdbview.addStyle({'within': {'distance': '5', 'sel': {'organic': True}}}, {'stick': {}})
    
    # Set surface for other atoms (except ligands) with VDW representation
    pdbview.addSurface(py3Dmol.VDW, {'opacity': 0.85, 'color': 'white'}, \
                       {'not': {'hetflag': True}})
    
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.zoom(0.8)  # Adjust zoom level as needed
    pdbview.spin(True, speed=0.025,)
    showmol(pdbview, height=500, width=800)


# Protein sequence input
DEFAULT_SEQ = "1sqt"
txt = st.sidebar.text_area('Input PDB-ID: with ONLY a lignad', DEFAULT_SEQ, height=275)


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


# PDB-fetching
def update(sequence=txt):
    pdb_id = sequence[:4]
    pdb_filename = get_pdb(pdb_id)

    with open(pdb_filename, 'r') as f:
        pdb_string = f.read()


    # Display protein structure
    st.subheader(f'Visualization of the {pdb_id} protein structure')
    render_mol(pdb_string)



predict = st.sidebar.button("Let's have a look!", on_click=update)


if not predict:
    st.warning('👈 Enter proten PDB-ID!')


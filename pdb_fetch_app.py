# This is app is created by Chanin Nantasenamat (Data Professor) https://youtube.com/dataprofessor
# Credit: This app is inspired by https://huggingface.co/spaces/osanseviero/esmfold

import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio

#st.set_page_config(layout = 'wide')
st.sidebar.title('🎈 Protein Preparation from RCSB.com')
st.sidebar.write('[*ESMFold*](https://esmatlas.com/about) is an end-to-end single sequence protein structure predictor based on the ESM-2 language model. For more information, read the [research article](https://www.biorxiv.org/content/10.1101/2022.07.20.500902v2) and the [news article](https://www.nature.com/articles/d41586-022-03539-1) published in *Nature*.')

# stmol
def render_mol(pdb):
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb,'pdb')
    pdbview.setStyle({'cartoon':{'color':'spectrum'}})
    pdbview.setBackgroundColor('white')#('0xeeeeee')
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    showmol(pdbview, height = 500,width=800)

# Protein sequence input
DEFAULT_SEQ = "6o0k"
txt = st.sidebar.text_area('Input PDB-ID', DEFAULT_SEQ, height=275)


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

    struct = bsio.load_structure(pdb_filename, extra_fields=["b_factor"])
    b_value = round(struct.b_factor.mean(), 4)

    # Display protein structure
    st.subheader('Visualization of predicted protein structure')
    render_mol(pdb_string)

    # plDDT value is stored in the B-factor field
    st.subheader('plDDT')
    st.write('plDDT is a per-residue estimate of the confidence in prediction on a scale from 0-100.')
    st.info(f'plDDT: {b_value}')

    st.download_button(
        label="Download PDB",
        data=pdb_string,
        file_name=f'{pdb_id}.pdb',
        mime='text/plain',
    )


predict = st.sidebar.button('Predict', on_click=update)


if not predict:
    st.warning('👈 Enter proten PDB-ID!')


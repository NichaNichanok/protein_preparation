# activate the virtual environment
#!source /home/nauevech/pp_env/bin/activate

import pandas as pd
import biopandas
from biopandas.pdb import PandasPdb
from prody import parsePDBHeader
from typing import Optional

def read_pdb_to_dataframe(pdb_path: Optional[str] = None,model_index: int = 1,parse_header: bool = True, ) -> pd.DataFrame:
    """
    Read a PDB file, and return a Pandas DataFrame containing the atomic coordinates and metadata.

    Args:
        pdb_path (str, optional): Path to a local PDB file to read. Defaults to None.
        model_index (int, optional): Index of the model to extract from the PDB file, in case
            it contains multiple models. Defaults to 1.
        parse_header (bool, optional): Whether to parse the PDB header and extract metadata.
            Defaults to True.

    Returns:
        pd.DataFrame: A DataFrame containing the atomic coordinates and metadata, with one row
            per atom
    """
    atomic_df = PandasPdb().read_pdb(pdb_path)
    if parse_header:
        header = parsePDBHeader(pdb_path)
    else:
        header = None
    atomic_df = atomic_df.get_model(model_index)
    if len(atomic_df.df["ATOM"]) == 0:
        raise ValueError(f"No model found for index: {model_index}")

    return pd.concat([atomic_df.df["ATOM"], atomic_df.df["HETATM"]]), header

pd = read_pdb_to_dataframe("/home/nauevech/Documents/protein_preparation/6o0k_stripped.pdb")
print(pd)

import requests
from bs4 import BeautifulSoup

def get_pdb(pdb_id):
    """
    Download the protein structure from the RCSB website by PDB ID
    
    Args:
    - PDB ID: 4 letter code protein from the RCSB.org website
    """
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(pdb_url)
    if response.status_code == 200:
        pdb_filename = f"{pdb_id}.pdb"
        with open(pdb_filename, 'wb') as f:
            f.write(response.content)
        return pdb_filename
    else:
        raise ValueError(f"Failed to download PDB file for {pdb_id}")


def fetch_ligand_name(pdb_id):
    """
    Search for the drug-ligand name defined in the PDB file 
    by fetching the information in the website, 
    i.e. the binding affinity annotations ID
    
    Args:
    - PDB ID: 4 letter code protein from the RCSB.org website
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


if __name__ == "__main__":
    # Example usage
    pdb_id = "1sqt"
    binding_id = fetch_ligand_name(pdb_id)
    if binding_id:
        print("Binding ID:", binding_id)
    else:
        print("Binding ID not found.")




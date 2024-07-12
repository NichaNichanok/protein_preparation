# Predict the Search space (grid box) for further molecular docking

- Translate the Af2bind colab notebook for local running: https://github.com/sokrypton/af2bind
- Based on binding residues predicted by the Af2bind
- Note: This test version (gridbox_af2bind.py) only worked for the single binding site at the moment






### Install aria2 and wget (for MacOS)
brew install aria2 wget

### Download AlphaFold parameters
aria2c -x 16 https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar
mkdir -p params
tar -xf alphafold_params_2021-07-14.tar -C params
touch params/done.txt

### Download af2bind parameters
wget https://github.com/sokrypton/af2bind/raw/main/attempt_7_2k_lam0-03.zip
unzip attempt_7_2k_lam0-03.zip -d af2bind_params

### Install Python packages
pip install numpy pandas jax scipy plotly py3Dmol matplotlib colabdesign

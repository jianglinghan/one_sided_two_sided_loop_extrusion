mamba create -n extrusion
mamba activate extrusion
mamba install -c conda-forge python=3.10.19
pip install openmm[cuda12]
pip install cuda-toolkit[all]==12.4
pip install -r requirements.txt
wget https://github.com/open2c/polychrom/archive/refs/tags/v0.1.0.tar.gz
tar xzvf v0.1.0.tar.gz
cd polychrom-0.1.0 && python setup.py install
wget https://github.com/open2c/polykit/archive/refs/tags/v0.0.0.tar.gz
tar xzvf v0.0.0.tar.gz
cd polykit-0.0.0 && python setup.py install 

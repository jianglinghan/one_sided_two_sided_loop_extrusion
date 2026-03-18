import os
import sys
import ast
import h5py
import glob
import joblib
import numpy as np
import cooler
from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file
import polychrom.contactmaps
from polykit.analysis import contact_maps as cms

indir, outdir = sys.argv[1:3]
# --- 配置参数 ---
RES = 2000  
N_BINS = 1000  
monomer_per_replica = 1000 
if not os.path.exists(outdir):
    os.makedirs(outdir)

def build_coolers():

    sim_folders = glob.glob(f"{indir}/folder*") 
    im_folders = sorted(sim_folders)
    
    for folder_path in sim_folders:
        sim_name = os.path.basename(folder_path).split('_')[-1]
        print(f"Processing {sim_name}...")
        
        try:
            URIs = polychrom.hdf5_format.list_URIs(folder_path)
            URIs_eq = [u for u in URIs if int(u.split("::")[-1]) >= 0]

            mrc = polychrom.contactmaps.monomerResolutionContactMapSubchains(
                URIs_eq,
                np.arange(0, 10000, monomer_per_replica),
                monomer_per_replica,
                cutoff=2.3, n=10
            )

            cool_uri = f'{outdir}/{sim_name}'
            cms.coolify(mrc, cool_uri, binsize=RES)
            print(f"Saved: {cool_uri}")
            #break
        except:
            pass

if __name__ == "__main__":
    build_coolers()

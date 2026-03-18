import os
import sys
import h5py 
import time
import joblib
import numpy as np
import pandas as pd

from lattice_translocators8 import LEFTranslocator, LEFTranslocatorDynamicBoundary
import funcs8 as funcs

import polychrom
from polychrom import polymerutils
from polychrom import forces
from polychrom import forcekits
from polychrom.simulation import Simulation
from polychrom.starting_conformations import grow_cubic
from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file
from polychrom.lib.extrusion import  bondUpdater
import warnings
import ast

filename = sys.argv[1]
paramfile = sys.argv[2]
out_dir = sys.argv[3]

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

params = [ast.literal_eval(i) for i in filename.split('folder_')[1].split('_')[1::2]]
face, back, clife, cof, life, slife, birth, pause, sep, site, \
    monomer, replica, steps, vel, _mode, switch, speed_mult, stability_factor, dwell_time = params

if int(_mode) == 1:
    mode = 'asymmetric'
elif int(_mode) == 2:
    mode = 'symmetric'

print(mode)

loop_id, CHROM, START, END, RES, CTCF_left_positions, CTCF_right_positions, site_types = joblib.load(paramfile)

paramdict={
            'CTCF_facestall':[face, face],
            'CTCF_backstall':[back, back],
            'CTCF_lifetime':[clife, clife],
            'CTCF_offtime':[cof, cof],
            'LEF_lifetime':[life, life],
            'LEF_stalled_lifetime':[slife, slife],
            'LEF_birth':[0.001, birth],
            'LEF_pause':[pause, pause],
            'LEF_separation':sep,
            'sites_per_monomer':site,
            'monomers_per_replica':monomer,
            'number_of_replica':replica,
            'steps':steps,
            'velocity_multiplier':vel
            }

paramdict_keys={
                'CTCF_facestall':'face',
                'CTCF_backstall':'back',
                'CTCF_lifetime':'Clife',
                'CTCF_offtime':'Cof',
                'LEF_lifetime':'life',
                'LEF_stalled_lifetime':'slife',
                'LEF_birth':'birth',
                'LEF_pause':'pause',
                'LEF_separation':'sep',
                'sites_per_monomer':'site',
                'monomers_per_replica':'monomer',
                'number_of_replica':'replica',
                'steps':'steps',
                'velocity_multiplier':'vel'
                }

file_name = funcs.paramdict_to_filename(paramdict, paramdict_keys)
folder = f'{out_dir}/' + 'folder_' + file_name.split('file_')[1] + f'_{loop_id}'
if os.path.exists(folder):
    print("already exist")
else:
    os.makedirs(folder)

monomers_per_replica = paramdict['monomers_per_replica'] 
sites_per_monomer = paramdict['sites_per_monomer']
sites_per_replica = monomers_per_replica * sites_per_monomer

########### 1d simulation parameters for lattice ###########
Trajn = 10000 
trajectory_length = Trajn * paramdict['sites_per_monomer'] 
num_dummy_steps = trajectory_length // 5 
blocksteps = 5 
bins = np.linspace(0, trajectory_length, blocksteps, dtype=int)
N = (paramdict['monomers_per_replica']*paramdict['number_of_replica'])
LEFNum = N // paramdict['LEF_separation']

# params for new model
paramdict['mode'] = mode
paramdict['switchProb'] = switch
paramdict['speedMultiplier'] = speed_mult
paramdict['stability_factor'] = stability_factor
paramdict['dwell_time'] = dwell_time

translocator = funcs.make_translocator(LEFTranslocatorDynamicBoundary,
                                 site_types,
                                 CTCF_left_positions,
                                 CTCF_right_positions, 
                                 **paramdict)

with h5py.File(folder+"/LEFPositions.h5", mode='w') as myfile:
    dset = myfile.create_dataset("positions", 
                                 shape=(trajectory_length, LEFNum, 2), 
                                 dtype=np.int32, 
                                 compression="gzip")
    
    translocator.steps(0)
    
    for st, end in zip(bins[:-1], bins[1:]):
        cur = []
        for i in range(st, end):
            if i % 1000 == 0:
                print(i)
            translocator.step()        
            cur.append(translocator.LEFs.copy())
        cur = np.array(cur)
        dset[st:end] = cur
    myfile.attrs["N"] = N * paramdict['sites_per_monomer']
    myfile.attrs["LEFNum"] = LEFNum

### Molecular dynamics simulaiton ###
myfile = h5py.File(folder + "/LEFPositions.h5", mode='r')
sites_per_monomer = paramdict['sites_per_monomer']
N = myfile.attrs["N"] // sites_per_monomer
print(N)
LEFNum = myfile.attrs["LEFNum"]
LEFpositions = myfile["positions"][::sites_per_monomer]// sites_per_monomer
Nframes = LEFpositions.shape[0]

stiff = 1
dens = 0.2
box = (N / dens) ** 0.33  

smcStepsPerBlock = 1  
data = grow_cubic(N, int(box) - 2)  
block = 0  
steps= paramdict['steps']


saveEveryBlocks = 10   
restartSimulationEveryBlocks = 100

smcBondWiggleDist = 0.2
smcBondDist = 0.5

assert (Nframes % restartSimulationEveryBlocks) == 0 
assert (restartSimulationEveryBlocks % saveEveryBlocks) == 0

savesPerSim = restartSimulationEveryBlocks // saveEveryBlocks
simInitsTotal  = (Nframes) // restartSimulationEveryBlocks 


tstp = 70 
tmst = 0.01

milker = polychrom.lib.extrusion.bondUpdater(LEFpositions)

reporter = HDF5Reporter(folder=folder, max_data_length=100, overwrite=True, blocks_only=False)

for iteration in range(simInitsTotal):
    a = Simulation(
            platform="cuda", GPU="0",
            integrator='langevin',  timestep=tstp, collision_rate=tmst,
            error_tol=0.01,  
            N = len(data),
            reporters=[reporter],
            PBCbox=[box, box, box],
            precision="mixed") 
    a.set_data(data)  

    a.add_force(
        forcekits.polymer_chains(
            a,
            chains=[(0, None, 0)],

            bond_force_func=forces.harmonic_bonds,
            bond_force_kwargs={
                'bondLength':1.0,
                'bondWiggleDistance':0.1, 
             },

            angle_force_func=forces.angle_force,
            angle_force_kwargs={
                'k':1.5
            },

            nonbonded_force_func=forces.polynomial_repulsive,
            nonbonded_force_kwargs={
                'trunc':1.5, 
                'radiusMult':1.05, 
            },
            except_bonds=True,
    ))
    kbond = a.kbondScalingFactor / (smcBondWiggleDist ** 2)
    bondDist = smcBondDist * a.length_scale

    activeParams = {"length":bondDist,"k":kbond}
    inactiveParams = {"length":bondDist, "k":0}
    milker.setParams(activeParams, inactiveParams)

    milker.setup(bondForce=a.force_dict['harmonic_bonds'],
                blocks=restartSimulationEveryBlocks)
    for t,l in enumerate(milker.allBonds):
        for b in l:
            if (b[0] == 11296) or (b[1] == 11296):
                print(t,b)
    if iteration==0:
        a.local_energy_minimization() 
    else:
        a._apply_forces()

    for i in range(restartSimulationEveryBlocks):        
        if i % saveEveryBlocks == (saveEveryBlocks - 1):  
            a.do_block(steps=steps)
        else:
            a.integrator.step(steps)  
        if i < restartSimulationEveryBlocks - 1: 
            curBonds, pastBonds = milker.step(a.context)  
    data = a.get_data()  
    del a

    reporter.blocks_only = True  

    time.sleep(0.2)  

reporter.dump_data()

myfile.close()









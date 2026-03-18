import sys
import joblib
import cooler
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm

data_dir = sys.argv[1]
paramfile = sys.argv[2] # cr139.pkl
title = 'test'
if len(sys.argv) > 3:
    title = sys.argv[3]

loop_id, CHROM, START, END, RES, CTCF_left_positions, CTCF_right_positions, site_types = joblib.load(paramfile)

start_bin = 300 
end_bin = 700 

loading_site = np.nonzero(site_types)[0]


cool_path = f'{data_dir}/{loop_id}.2000.cool'

ctcf_left_plot = sorted(list(set(CTCF_left_positions // 10))) 
ctcf_right_plot = sorted(list(set(CTCF_right_positions // 10)))
loading_plot = sorted(list(set(loading_site // 10)))

c = cooler.Cooler(cool_path)
mat = c.matrix(balance=False, sparse=False).fetch(c.chromnames[0])
mat = mat.astype(float)
mat /= np.median(np.diag(mat, 2))

fig, axes = plt.subplots(2, 1, gridspec_kw={
                             'height_ratios': [5, 2], 
                             'hspace': 0.1
                         })

ax_sim = axes[0]
ax_track = axes[1]

mat = mat[start_bin:end_bin, start_bin:end_bin]
vmax = np.percentile(mat, 97)
vmin = np.percentile(mat, 5)
vmax, vmin = 0.0521552931759054, 0.005
CMAP = LinearSegmentedColormap.from_list('interaction', ['#FFFFFF','#FFDFDF','#F70000'])
im = ax_sim.imshow(mat, norm=LogNorm(vmin, vmax), cmap=CMAP)

ax_sim.set_title(f"Simulation: {title}")

cbar = plt.colorbar(im, ax=ax_sim, fraction=0.046, pad=0.04, shrink=0.3)        

ax_track.scatter(ctcf_left_plot, [1]*len(ctcf_left_plot), s=10, marker='>', color='green', label='CTCF R')
ax_track.scatter(ctcf_right_plot, [0]*len(ctcf_right_plot), s=10, marker='<', color='blue', label='CTCF L')
ax_track.scatter(loading_plot, [2]*len(loading_plot), s=10, c='red', marker='|', label='Loading')

ax_track.set_yticks([0, 1, 2])
ax_track.set_yticklabels(['CTCF R', 'CTCF L', 'Loading'], fontsize=5)
ax_track.set_xlim(start_bin, end_bin) 
ax_track.set_ylim(-1, 3)
ax_track.grid(axis='x', linestyle='--', alpha=0.3)
pos = ax_sim.get_position()
ax_track.set_position([pos.x0, pos.y0 - 0.1, pos.width, 0.05])

plt.savefig(f'{title}.{loop_id}_sim_only.zoom.svg', bbox_inches='tight')
print(f"Plot saved as {title}.{loop_id}_sim_only.zoom.svg")

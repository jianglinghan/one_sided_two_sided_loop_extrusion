# Simulation Framework for One-sided and Two-sided Loop Extrusion Models [![DOI](https://zenodo.org/badge/1185208969.svg)](https://doi.org/10.5281/zenodo.19093828)
<img width="1313" height="490" alt="title" src="https://github.com/user-attachments/assets/506a7975-4296-4a22-b909-2eac47829db2" />


## Description
This repository contains the source code and configuration files for simulating distinct loop extrusion models. It is designed to characterize the dynamics of chromatin loop formation and maintenance under different mechanistic assumptions.

Specifically, we implemented two models. Accordingly, the repository is organized into two primary folders:
* **`two-sided`**: Implementing the classical two-sided extrusion model.
* **`one-sided`**: Implementing the one-sided, two-phase loop extrusion model.

Each folder contains the following core components for simulation and post-processing:

### 1. Simulation Core
* **`simu_all.py`**: The main simulation engine.
* **`lattice_translocators8.py`**: Implements the dynamics of loop-extruding factors (LEFs) on the 1D lattice.
* **`funcs8.py`**: Helper functions used by `lattice_translocators8.py`.
* **`cr139.pkl`**: Example input dataset containing the following information in order: `loop_id`, `CHROM`, `START`, `END`, `RES`, `CTCF_left_positions`, `CTCF_right_positions`, `LEF_sites`.

### 2. Execution Scripts
* **`step1-simu.sh`**: Runs the loop extrusion simulation using `simu_all.py`.
* **`step2-build_cools.sh`**: Converts simulation trajectories into chromatin contact maps in `.cool` format using `build_cools.py`.
* **`step3-eg_simu.sh`**: Generates the whole-region and the zoom-in heatmaps for the demo dataset using `eg_simu.py` and `eg_zoom.py`.

---

## Requirements
* **Python** (tested on 3.10.19)
* **[polychrom](https://github.com/open2c/polychrom)**: toolkit for polymer simulations.
* **[openMM](https://github.com/openmm/openmm)**: molecular dynamics engine. Please make sure CUDA and OpenMM versions are compatible when using NVIDIA GPUs. Refer to [OpenMM documentation](https://docs.openmm.org/latest/userguide/application/01_getting_started.html) for more details.
* **[polykit](https://github.com/open2c/polykit)**: toolkit for polymer simulation setup and analysis.
* Other required packages and versions are listed in `requirements.txt`.

---
## Installation
Users can create a simulation environment using Conda or Mamba by following these steps in sequence:
```bash
mamba create -n extrusion
mamba activate extrusion
mamba install -c conda-forge python=3.10.19
# If an NVIDIA GPU is available:
# pip install openmm[cuda12]
# pip install cuda-toolkit[all]==12.4
# else
# pip install openmm
pip install -r requirements.txt
wget https://github.com/open2c/polychrom/archive/refs/tags/v0.1.0.tar.gz
tar xzvf v0.1.0.tar.gz
cd polychrom-0.1.0 && python setup.py install
wget https://github.com/open2c/polykit/archive/refs/tags/v0.0.0.tar.gz
tar xzvf v0.0.0.tar.gz
cd polykit-0.0.0 && python setup.py install
```
---

## Workflow

1.  **Enter a model directory.** For example, to run the one-sided, two-phase loop extrusion model:
    ```bash
    cd one-sided
    ```

2.  **Run simulation**
    ```bash
    python simu_all.py <simulation paramters> <LEF_CTCF_position_file> <simulation_output_dir>
    ```
    **Example:**
    ```bash
    python simu_all.py folder_face_1_back_0_Clife_100_Cof_300_life_200_slife_1200_birth_0.05_pause_0.5_sep_50_site_10_monomer_1000_replica_10_steps_155_vel_1_mode_1_switch_1_speed_2_stability_0.01_dwell_200 cr139.pkl wild
    ```
    > If an NVIDIA GPU is unavailable, you can use `simu_all.cpu.py` instead of `simu_all.py`, though it will be less efficient.

3.  **Generate contact maps (.cool):**
    ```bash
    python build_cools.py <simulation_output_dir> <cool_output_dir>
    ```
    **Example:**
    ```bash
    python build_cools.py wild wild/cool_outputs
    ```

4.  **Visualization**

    **Whole-region heatmap:**
    ```bash
    python eg_simu.py <cool_output_dir> <LEF_CTCF_position_file> <title>
    ```
    **Example:**
    ```bash
    python eg_simu.py wild/cool_outputs cr139.pkl wild
    ```

    **Zoom-in heatmap:**
    ```bash
    python eg_zoom.py <cool_output_dir> <LEF_CTCF_position_file> <title>
    ```
    **Example:**
    ```bash
    python eg_zoom.py wild/cool_outputs cr139.pkl wild
    ```


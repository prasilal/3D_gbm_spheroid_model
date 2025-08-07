
# 3D GBM Spheroid Model

This repository contains a biologically motivated 3D spheroid model of glioblastoma (GBM) spheroid growth under co-culture conditions with mesenchymal stromal cells (MSCs). The framework is based on a **Cellular Potts Model (CPM)** with integrated **reaction–diffusion dynamics**, **logic-based cell state transitions**, and **interactive visualization**.

## Model Summary

- 3D Cellular Potts Model
- Multiple cell types: GBM, glioblastoma stem-like cells (GSC), and MSCs
- Cell state transitions: Proliferation, quiescence, necrosis, de/differentiation
- Field interactions: Glucose, oxygen, lactate, ECM, growth signals
- Michaelis-Menten nutrient uptake
- Custom diffusion and decay per field and cell type
- GLMakie 3D visualizations of cells, fields, and simulation statistics

## Repository Structure

```
root/
├── analysis/                           # Jupyter notebooks for (sim and exp) data analysis
│   ├── cell_size_analysis.ipynb
│   ├── intensity_analysis.ipynb
│   └── ... (other analysis notebooks)
│
├── data/                               # Experimental and simulation data
│   ├── experimental_data/                     # Raw experimental data
│   ├── experimental_data_analysis/            # Processed analysis results from experimental data
│   ├── processed/                             # Processed stats from simulation data analysis
│   └── raw_simulation_outputs/                # .jld2 files from simulations
│
├── media/                              # All output media
│   ├── thesis_figures/                        # Figures used in the thesis
│   └── other/                                 # Supplementary plots or visuals
│
├── model/                              # Simulation codebase (Julia)
│   ├── main/                                  # Main simulation code (used for development and plotting)
│   │   ├── scripts/                                  # Entry point scripts for different simulation runs
│   │   │   ├── 3D_h-model.jl                                # Main model file
│   │   │   ├── example_cfg.jl                               # Example cfg with description and tutorial comments
│   │   │   ├── ... (other configuration scripts)
│   │   │   ├── hooks/                                       # Additional logic hooks (division, state transition, outer actions, stats and visualization,...)
│   │   │   └── tests/                                       # Tests or small-scale examples
│   │   ├── src/                                      # Underlying simulation logic 
│   │   ├── Project.toml
│   │   └── Manifest.toml
│   │
│   └── headless/                              # Headless (non-visual) version for cluster runs
│       ├── scripts/                                  # Batch execution or cluster-friendly variants
│       ├── src/ 
│       ├── Project.toml
│       └── Manifest.toml

```

## Installation

### Julia Setup

Requires **Julia ≥ 1.9.0**

1. Clone the repository and navigate to the `model/main/scripts` folder:

    ```bash
    git clone https://github.com/prasilal/gbm-spheroid-model.git
    cd model/main/scripts
    ```

2. Start Julia and activate the environment:

    ```julia
    using Pkg
    Pkg.activate("..")
    Pkg.instantiate()
    ```

If you plan to run simulations on a cluster, repeat these steps in `model/headless/`.

### Python Setup (for analysis)

```bash
pip install -r requirements.txt
```

## Running a Simulation

### Main Configurable Simulation

To run the primary model with customizable settings, use `test_run.jl`, which internally loads the full simulation logic from `3D_h-model.jl` and runs a parameterized version of `example_cfg.jl`:

```julia
include("test_run.jl")
```

You can customize parameters like grid size, cell ratios, initial seed, or number of simulation steps directly inside `test_run.jl`.
This is the recommended entry point for running full-scale or thesis-relevant simulations.

### Running Test or Small-Scale Simulations

Smaller test runs, visualizations, or experimental logic (e.g., division tests, field visualization, parameter tuning) can be run by first including the `3D_h-model.jl` and then simply including the appropriate script:

```julia
include("absorb_test.jl")  
```

These are typically located in:

```
model/main/scripts/tests/
```

or alongside the main config scripts in:

```
model/main/scripts/
```

### Cluster Execution

To run on a cluster, use the scripts in:

```
model/headless/scripts/
```

These are designed to run simulations in batch without GUI and to export `.jld2` output at customizable intervals. You can modify these scripts to fit batch jobs or parameter sweeps (see `job.sh` for an example).

## Output

Simulation results can include:

- `.jld2` files containing CPM cell states, voxel field data, and grid configurations
- Rendered plots of 3D spheroid structure, nutrient fields, and cell state dynamics
- Post-analysis output from Python scripts

## Visualization

GLMakie is used for rendering:

- 3D volumetric views of the spheroid
- Sliced field visualizations (e.g. oxygen, lactate)
- Overlay statistics panel with live cell counts, state distributions, etc.

All visualization code is located in `model/main/scripts/3D_h-model.jl`.

## Documentation

The simulation code is fully documented with Julia docstrings. Major files:

| File | Description |
|------|-------------|
| `3D_h-model.jl` | Main simulation loop, config parsing, simulation logic |
| `transition_hooks.jl` | Modular state transition dispatcher (per cell type) |
| `division_hooks.jl` | Division and cell cycle logic |
| `example_cfg.jl` | Example configuration file with documented logical modules |

## Acknowledgments

This work was developed as part of a master’s thesis project on 3D GBM spheroid modeling. Simulations were run locally and on the **MetaCentrum** computing cluster, supported by the e-INFRA CZ project (ID:90254), funded by the Ministry of Education, Youth and Sports of the Czech Republic.

## License

MIT License — Free to use, modify, and redistribute with attribution.


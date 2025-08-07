# test_run.jl — minimal script to run a simulation with visualization on a local machine

# Define simulation arguments (can be edited manually)
# - "simscript": path to the file that defines mk_cfg(args) and run_sim(args)
# - "seed": random seed for reproducibility
# - "gridsize": cubic simulation grid size (X × Y × Z)
# - "nofcells": total number of cells to initialize
# - "nofsteps": number of macrosteps to run
# - "gbm", "gsc", "msc": fractions of each cell type (must sum to 1.0)


args = Dict(
    "simscript" => "example_cfg.jl",
    "seed" => 123,
    "gridsize" => 70,
    "nofcells" => 1000,
    "nofsteps" => 10,
    "gbm" => 0.735,
    "gsc" => 0.015,
    "msc" => 0.25
)

# Include the main simulation logic
include("3D_h-model.jl")

# Include the selected simulation configuration file
# This file must define mk_cfg(args) and run_sim(args)
simscript = args["simscript"]
include(simscript)

# Run the simulation using the given arguments
run_sim(args)
include("plots_stats_vis_hooks.jl")
include("division_hooks.jl")
include("transition_hooks.jl")
include("outer_actions_hooks.jl")


using JLD2

function save_sims_state(sim, step::Int)

    if sim.model.state.t <= sim.cfg[:runtime][:burnin_steps]
        return
    end
    folder = "sim_output"
    filename = "step_" * lpad(string(step), 5, '0') * ".jld2"
    filepath = joinpath(folder, filename)

    mkpath(folder)

    jldopen(filepath, "w") do f
        f["step"] = step

        # Save cell_voxels as group of arrays
        for (id, vox) in get_stat(sim.model, :pixels_by_cell)
            id > 2 || continue
            f["cell_voxels/$(id)"] = collect(Tuple.(vox))  # Nx3 array
        end

        # Save cell_states as group of strings
        for c in sim.cells
            if c.id > 2
                f["cell_states/$(c.id)"] = string(convert_idx_to_state(c.state[:cum_state], sim))
            end
        end
    end
    println("step $step saved to $filepath, but timestep: $(sim.model.state.t)")
    
    if step % 10 == 0
        save_vaxes_state(sim, step)
    end
end

function save_vaxes_state(sim, step::Int)
    folder = joinpath("sim_output", "vaxes")
    mkpath(folder)

    filename = "vaxes_step_" * lpad(string(step), 5, '0') * ".jld2"
    filepath = joinpath(folder, filename)

    jldopen(filepath, "w") do f
        f["step"] = step
        
        for (vname, vcfg) in sim.cfg[:vaxes]
            vid = vcfg[:id]
            grid = sim.vaxes[vid]

            f["vaxes/$(vname)"] = Array(grid)
        end
    end

    println("vaxes saved to $filepath")
end


function holder(sim, step::Int)
    # Extract cell voxel mappings, excluding barrier cells (id ≤ 2)
    cell_voxels = Dict(k => v for (k, v) in get_stat(sim.model, :pixels_by_cell) if k > 2)

    # Extract cell states, converted to string, also excluding barrier cells
    cell_states = Dict(c.id => string(convert_idx_to_state(c.state[:cum_state], sim))
                        for c in sim.cells if c.id > 2)

    # Build filename and path
    filename = "step_" * lpad(string(step), 5, '0') * ".jld2"
    filepath = joinpath("sim_output", filename)

    # Save the data
    @save filepath cell_voxels cell_states step
end

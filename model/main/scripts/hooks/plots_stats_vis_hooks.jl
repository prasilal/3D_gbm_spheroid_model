"""
plots_stats_vis_hook.jl

Provides visualization and analysis utilities for simulation results.

Includes:
- Cell state tracking over time (line plots, Sankey diagram)
- Field profile and slice plotting (XY slices, X-line, depth profiles)
- Utility functions for saving visualizations

Used in:
- Simulation visualization
- Debugging and analysis
- Figure generation for publications
"""

# Global tracker for cell state counts per step
const state_counts_log = []

# Internal counter to keep track of simulation steps logged
const step_counter = Ref(0)


"""
log_cell_states(sim, cell, new_state_idx, weight, _)

Hook function for tracking the number of cells in each state over time.

Triggered on state transition events, but logs only once per step
(by checking for `cell.id == 3`).

Adds a dictionary of state counts to `state_counts_log`.
"""
function log_cell_states(sim, cell, new_state_idx, weight, _)

    if cell.id == 3
        counts = Dict{Symbol, Int}()
        for c in sim.cells
            state_sym = convert_idx_to_state(c.state[:cum_state], sim)
            counts[state_sym] = get(counts, state_sym, 0) + 1
        end
        push!(state_counts_log, counts)
        step_counter[] += 1
    end
end


"""
plot_state_counts(save_path="outputs/state_counts.png")

Plots the temporal dynamics of cell state counts using data in `state_counts_log`.

Each state is shown as a separate line; counts are collected per simulation step.
"""
function plot_state_counts(save_path::String = "outputs/state_counts.png")
    states = unique(reduce(vcat, collect.(keys.(state_counts_log))))
    timepoints = 1:length(state_counts_log)

    fig = Figure(size = (800, 600)) 
    ax = GLMakie.Axis(fig[1, 1], xlabel = "Time step", ylabel = "Cell count", title = "Cell state dynamics")

    for state in sort(states)
        counts = [get(frame, state, 0) for frame in state_counts_log]
        lines!(ax, timepoints, counts, label = String(state))
    end

    axislegend(ax)
    save(save_path, fig)
end

"""
prepare_multistep_sankey_data(transition_log)

Converts a log of state transitions into node/edge data
for visualizing temporal transitions as a Sankey diagram.

Returns:
- `node_labels`: All states + steps as strings
- `links`: Map of transition tuples to flow counts
"""
function prepare_multistep_sankey_data(transition_log)
    node_labels = Set{String}()
    links = Dict{Tuple{String, String}, Int}()

    for trans in transition_log
        # Source = from state at current step
        src = "$(trans.from) at step $(trans.step)"
        # Target = to state at next step
        tgt = "$(trans.to) at step $(trans.step + 1)"

        # Register nodes
        push!(node_labels, src)
        push!(node_labels, tgt)

        # Count flow
        key = (src, tgt)
        links[key] = get(links, key, 0) + 1
    end

    return collect(node_labels), links
end


"""
parse_node_label(label)

Parses labels of the form "state at step X" and returns:
- state::String
- step::Int
"""
function parse_node_label(label::String)
    m = match(r"(.*) at step (\d+)", label)
    if m === nothing
        error("Could not parse label: $label")
    end
    return m.captures[1], parse(Int, m.captures[2])
end


"""
plot_multistep_sankey()

Draws a Sankey-like diagram from transition log,
showing how cells change states across steps.
Width of line encodes number of transitions.
"""
function plot_multistep_sankey()
    node_labels, links = prepare_multistep_sankey_data(transition_log)

    parsed_nodes = [parse_node_label(label) for label in node_labels]

    states = sort(unique(s[1] for s in parsed_nodes))
    state_to_y = Dict(state => 1.0 - (i-1)/(length(states)-1) for (i, state) in enumerate(states))

    steps = sort(unique(s[2] for s in parsed_nodes))
    max_step = maximum(steps)

    fig = Figure(size = (1400, 800))
    ax = GLMakie.Axis(fig[1, 1],
        xlabel = "Step",
        ylabel = "State",
        xticks = (steps, string.(steps)),
        yticks = (collect(values(state_to_y)), states),
        xgridvisible = true,
        ygridvisible = false
    )

    hidespines!(ax, :t, :r, :b)

    # Group transitions by step for better label offsets
    transitions_by_step = Dict{Int, Vector{Tuple{Tuple{String, String}, Int}}}()
    for (link, value) in links
        src, tgt = link
        src_state, src_step = parse_node_label(src)
        push!(get!(transitions_by_step, src_step, Vector{Tuple{Tuple{String, String}, Int}}()), (link, value))
    end

    # Draw transitions
    for (step, transitions) in transitions_by_step
        n = length(transitions)
        for (i, ((src_label, tgt_label), value)) in enumerate(transitions)
            src_state, src_step = parse_node_label(src_label)
            tgt_state, tgt_step = parse_node_label(tgt_label)

            x0 = src_step
            x1 = tgt_step
            y0 = state_to_y[src_state]
            y1 = state_to_y[tgt_state]

            lines!(ax, [x0, x1], [y0, y1]; linewidth = value * 0.4, color = :blue, transparency = true)

            mid_x = (x0 + x1) / 2
            mid_y = (y0 + y1) / 2

            x_offset = (i - (n+1)/2) * 0.05

            text!(ax, mid_x - x_offset, mid_y + 0.02, text = string(value), align = (:center, :bottom), fontsize=15)
        end
    end

    display(fig)
end


"""
print_average_volumes_by_type(sim)

Prints average voxel volume per cell type based on CPM cell occupancy.

Useful for debugging or validating cell growth and volume constraints.
"""
function print_average_volumes_by_type(sim)
    volumes = Dict{Symbol, Vector{Int}}()

    for cell in values(sim.cells)
        ctype = cell.state[:type]
        vol = length(get(get_stat(sim.model, :pixels_by_cell), cell.id, []))
        print("$ctype: ")
        println(vol)
        push!(get!(volumes, ctype, []), vol)
    end

    println("Average cell volumes by type:")
    for (ctype, vols) in volumes
        avg = mean(vols)
        println("  $ctype: $(round(avg, digits=2)) voxels")
    end
end

"""
plot_outflow_history(vax::Symbol, save_path="outflow.png")

Plots the time series of average outflow (absorption) for a given field (`:glucose` or `:oxygen`).

Arguments:
- `vax`: Field name (`:glucose` or `:oxygen`)
- `save_path`: Optional output path for saving the plot (default "outflow.png")

Uses:
- `glucose_outflow_tracker` or `oxygen_outflow_tracker`
"""
function plot_outflow_history(vax::Symbol, save_path::String = "outflow.png")
    if vax == :glucose
        history = glucose_outflow_tracker
        title = "Glucose Outflow History"
    elseif vax == :oxygen
        history = oxygen_outflow_tracker
        title = "Oxygen Outflow History"
    else
        error("Unknown vax type: $vax")
    end

    fig = Figure(size = (800, 400))
    ax = GLMakie.Axis(fig[1, 1]; title=title, xlabel="Time step", ylabel="Outflow")
    lines!(ax, 1:length(history), history)
    display(fig)
    save(save_path, fig)
end

"""
plot_averaged_x_profile(field::Array{Float64,3})

Plots the average field value along X, averaged across YZ for each X slice.
Useful for showing depth gradients or nutrient penetration profiles.
"""
function plot_averaged_x_profile(field::Array{Float64,3})
    x_vals = 1:size(field, 1)
    profile = [mean(field[x, :, :]) for x in x_vals]
    println("im here")

    fig = Figure(size = (800, 400))
    ax = GLMakie.Axis(fig[1, 1]; xlabel="x", ylabel="Average Field", title="Mean across YZ plane")
    lines!(ax, x_vals, profile)
    display(fig)
end

"""
plot_line_profile(field::Array{Float64,3}, y::Int, z::Int)

Plots the field value along the X-axis at fixed Y and Z coordinates.
Allows examining a specific cross-section of the field in detail.
"""
function plot_line_profile(field::Array{Float64,3}, y::Int, z::Int)
    x_vals = 1:size(field, 1)
    profile = [field[x, y, z] for x in x_vals]

    fig = Figure(size = (800, 400))
    ax = GLMakie.Axis(fig[1, 1]; xlabel="x", ylabel="Field value", title="Line profile at y=$y, z=$z")
    lines!(ax, x_vals, profile)
    display(fig)
end

"""
plot_field_slice_xy(field, z_plane::Int)

Displays a heatmap of the field in the XY plane at depth `z_plane`.
Uses a viridis colormap and shows raw field values.
"""
function plot_field_slice_xy(field, z_plane::Int)
    fig = Figure(size = (600, 500))
    ax = GLMakie.Axis(fig[1, 1], title = "Field @ z = $z_plane", aspect = DataAspect())

    data = field[:, :, z_plane]'
    println(data)
    hm = heatmap!(ax, data, colormap = :viridis)

    Colorbar(fig[1, 2], hm)
    display(fig)
end

"""
save_field_slice_xy(field, z_plane::Int; field_name="Field", save_path="field_slice.png")

Saves a heatmap image of the XY plane at `z_plane` for the given 3D field.

Arguments:
- `field_name`: Used in the title and colorbar label
- `save_path`: Output file name (default: "field_slice.png")
"""
function save_field_slice_xy(field, z_plane::Int; field_name::String = "Field", save_path::String = "field_slice.png")
    fig = Figure(size = (700, 550))
    ax = GLMakie.Axis(fig[1, 1],
        title = "$field_name concentration at z = $z_plane",
        xlabel = "X (voxels)",
        ylabel = "Y (voxels)",
        aspect = DataAspect()
    )

    data = field[:, :, z_plane]'
    hm = heatmap!(ax, data, colormap = :viridis)
    Colorbar(fig[1, 2], hm, label = "$field_name")

    save(save_path, fig)
end

"""
save_spheroid_slice_image(sim, z=50, filename="spheroid_slice.png")

Saves a colored slice of the cell population at a given Z plane (`z`).
Each cell type is assigned a numerical color:
- GBM: red
- GSC: blue
- MSC: green

The function uses per-voxel information from the CPM to color each pixel.
"""
function save_spheroid_slice_image(sim, z=50, filename="spheroid_slice.png")
    grid_size = size(sim.model.grid)
    slice_img = fill(0, grid_size[1], grid_size[2])  # slice at z

    # Assign numerical color codes for each type
    type_colors = Dict(:gbm => 1, :gsc => 2, :msc => 3)

    for cell in values(sim.cells)
        ctype = cell.state[:type]
        color_val = get(type_colors, ctype, 0)
        cell_pixels = get(get_stat(sim.model, :pixels_by_cell), cell.id, [])
        for idx in cell_pixels
            x, y, zz = Tuple(idx)
            if zz == z
                slice_img[x, y] = color_val
            end
        end
    end

    fig = Figure(size=(600, 600))
    ax = GLMakie.Axis(fig[1, 1]; title="Spheroid Slice at z = $z", aspect=DataAspect())
    heatmap!(ax, slice_img'; colormap = [:white, :red, :blue, :green], interpolate=false)
    save(filename, fig)
    println("Slice image saved to $filename")
end
"""
outer_actions_hooks.jl

Defines external interventions in the simulation grid (e.g., injection and absorption of fields)
and utilities for visualizing chemical field dynamics over time and space.

Main Features:
- Inject or absorb oxygen and glucose at boundary regions
- Track and plot nutrient outflow trends
- Visualize field slices or profiles across spatial dimensions
- Export spheroid structure as z-slices
"""


const oxygen_outflow_tracker = []
const glucose_outflow_tracker = []


"""
inject_absorb_track_function(cfg, model, r, cid, center)

Applies an outer action (injection or absorption) over a region `r` during the simulation.

Arguments:
- `cfg`: Configuration dictionary (must contain `:sim`, `:outer_action`, etc.)
- `model`: GridModel object from the CPM
- `r`: Region to act on (defined as ((x1,y1,z1), (x2,y2,z2), :inject|:absorb))
- `cid`, `center`: Unused in this function but required for compatibility with outer action API

Behavior:
- If `r[3] == :inject`: sets glucose/oxygen levels in the region to predefined injection values
- If `r[3] == :absorb`: zeros the field and tracks the removed amount per timestep

Tracks:
- `oxygen_outflow_tracker` and `glucose_outflow_tracker` record per-step average outflow

Note:
- Does nothing during the burn-in period
"""
function inject_absorb_track_function(cfg, model, r, cid, center)
    sim = cfg[:sim]
    # Skip during burn-in
    if model.state.t <= cfg[:runtime][:burnin_steps]
        return
    end
    
    # Prepare access to chemical fields
    # glucose
    glucose_field = sim.vaxes[cfg[:vaxes][:glucose][:id]].x

    # oxygen
    oxygen_field = sim.vaxes[cfg[:vaxes][:oxygen][:id]].x
    x_range = r[1][1]:r[2][1]
    y_range = r[1][2]:r[2][2]
    z_range = r[1][3]:r[2][3]

    total_out_glucose = 0.0
    total_out_oxygen = 0.0
    if r[3] == :inject
    end
    if r[3] == :absorb
    end

    for x in x_range, y in y_range, z in z_range
        if r[3] == :inject
            glucose_field[CartesianIndex(x,y,z)] = cfg[:outer_action][:injection_amount][:glucose]
            oxygen_field[CartesianIndex(x,y,z)] = cfg[:outer_action][:injection_amount][:oxygen]
        end
        if r[3] == :absorb
            total_out_glucose += glucose_field[CartesianIndex(x,y,z)]
            glucose_field[CartesianIndex(x,y,z)] = 0.0
            total_out_oxygen += oxygen_field[CartesianIndex(x,y,z)]
            oxygen_field[CartesianIndex(x,y,z)] = 0.0
        end
    end
    if r[3] == :absorb 
        average_out_glucose = total_out_glucose / length(x_range) / length(y_range) / length(z_range)
        average_out_oxygen = total_out_oxygen / length(x_range) / length(y_range) / length(z_range)
        push!(glucose_outflow_tracker, average_out_glucose)
        push!(oxygen_outflow_tracker, average_out_oxygen)
    end

end



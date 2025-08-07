
const oxygen_outflow_tracker = []
const glucose_outflow_tracker = []

function inject_absorb_track_function(cfg, model, r, cid, center)
    sim = cfg[:sim]
    if model.state.t <= cfg[:runtime][:burnin_steps]
        return
    end
    
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
    #println("Timestep: $(model.state.t- cfg[:runtime][:burnin_steps])")

end



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


function plot_averaged_x_profile(field::Array{Float64,3})
    x_vals = 1:size(field, 1)
    profile = [mean(field[x, :, :]) for x in x_vals]
    println("im here")

    fig = Figure(size = (800, 400))
    ax = GLMakie.Axis(fig[1, 1]; xlabel="x", ylabel="Average Field", title="Mean across YZ plane")
    lines!(ax, x_vals, profile)
    display(fig)
end

function plot_line_profile(field::Array{Float64,3}, y::Int, z::Int)
    x_vals = 1:size(field, 1)
    profile = [field[x, y, z] for x in x_vals]

    fig = Figure(size = (800, 400))
    ax = GLMakie.Axis(fig[1, 1]; xlabel="x", ylabel="Field value", title="Line profile at y=$y, z=$z")
    lines!(ax, x_vals, profile)
    display(fig)
end

function plot_field_slice_xy(field, z_plane::Int)
    fig = Figure(size = (600, 500))
    ax = GLMakie.Axis(fig[1, 1], title = "Field @ z = $z_plane", aspect = DataAspect())

    data = field[:, :, z_plane]'
    println(data)
    hm = heatmap!(ax, data, colormap = :viridis)

    Colorbar(fig[1, 2], hm)
    display(fig)
end


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
    println("✅ Slice image saved to $filename")
end


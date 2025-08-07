module Visualization

const RGBAf0 = RGBA{Float32}

function hex2rgba(hex)
    RGBA(((hex & 0xff0000) >> 16) / 255,
         ((hex & 0xff00) >> 8) / 255,
         (hex & 0xff) / 255,
         1.0)
end

rc(hex) = ((hex & 0xff0000) >> 16) / 255
gc(hex) = ((hex & 0xff00) >> 8) / 255
bc(hex) = (hex & 0xff) / 255

default_colors = Dict( # default colors for cells
    0 => hex2rgba(0x370617), # cell border
    1 => RGBA(1,1,1,1), # bg
    2 => hex2rgba(0x9fa6a6),
    3 => hex2rgba(0xf1f3c2),
    4 => hex2rgba(0xf3e700),
    5 => hex2rgba(0xf3d368),
    6 => hex2rgba(0xf3ba9e),
    7 => hex2rgba(0xf37b35),
    8 => hex2rgba(0xf31000),
    9 => hex2rgba(0xf3b4f1),
    10 => hex2rgba(0xf300ae),
    11 => hex2rgba(0x9a00f3),
    12 => hex2rgba(0x6100f3),
    13 => hex2rgba(0xb48bf3),
    14 => hex2rgba(0x648ff3),
    15 => hex2rgba(0xadd6f3),
    16 => hex2rgba(0x1119f3),
    17 => hex2rgba(0x1df3e8),
    18 => hex2rgba(0x20f386),
    19 => hex2rgba(0xbff3a3),
    20 => hex2rgba(0x00a876),
    21 => hex2rgba(0x961200),
    22 => hex2rgba(0xffeedb)
)

default_colors_field = Dict( # default colors for vaxes
    1 => hex2rgba(0xf1f3c2),
    2 => hex2rgba(0xf3e700),
    3 => hex2rgba(0xf3d368),
    4 => hex2rgba(0xf3ba9e),
    5 => hex2rgba(0xf37b35),
    6 => hex2rgba(0xf31000),
    7 => hex2rgba(0xf3b4f1),
    8 => hex2rgba(0xf300ae),
    9 => hex2rgba(0x9a00f3),
    10 => hex2rgba(0x6100f3),
    11 => hex2rgba(0xb48bf3),
    12 => hex2rgba(0x648ff3),
    13 => hex2rgba(0xadd6f3),
    14 => hex2rgba(0x1119f3),
    15 => hex2rgba(0x1df3e8),
    16 => hex2rgba(0x20f386),
    17 => hex2rgba(0xbff3a3),
    18 => hex2rgba(0x00a876),
    19 => hex2rgba(0x961200),
    20 => hex2rgba(0xffeedb)
)

function draw_cpm(sim;
                  cell_kind_color = default_colors,
                  num_cell_kinds = nothing,
                  show_activity = false,
                  show_fields = false,
                  field_color = default_colors_field,
                  vox = nothing)
    cpm = sim.model
    if show_fields
        for v in values(sim.cfg[:vaxes])
            if get(v,:show,false)
                # vox = draw_field_contour(sim.vaxes[v[:id]],
                #                          haskey(v,:color) ? v[:color] : field_color[v[:id]],
                #                          vox = vox, nsteps = 30)
                vox = draw_field(sim.vaxes[v[:id]],
                                 haskey(v,:color) ? v[:color] : field_color[v[:id]],
                                 vox = vox)
            end
        end
    end

    vox = draw_cells(cpm, cid -> cell_kind_color[cell_kind(cpm, cid)], vox = vox)

    num_cell_kinds = num_cell_kinds == nothing ? maximum(cpm.state.cell_to_kind) : num_cell_kinds

    for k in 2:num_cell_kinds
        draw_cell_borders(cpm, k, cell_kind_color[0], vox = vox)
        if show_activity
            draw_activity_values(cpm, k, vox = vox)
        end
    end

    vox
end

function draw_cpm_3d(sim;
                     cell_kind_color = default_colors,
                     show_fields = true,
                     show_cells = true,
                     show_activity = false,
                     field_color = default_colors_field)

    fig = Figure(size = (1200, 600))
    ax1 = Axis3(fig[1, 1], title = "Cells", aspect = :data)
    ax2 = Axis3(fig[1, 2], title = "Vaxes", aspect = :data)

    if show_cells
        draw_cells_3d!(ax1, sim, cell_kind_color)
    end

    if show_fields
        for v in values(sim.cfg[:vaxes])
            if get(v, :show, false)
                field = sim.vaxes[v[:id]]
                col = get(v, :color, field_color[v[:id]])
                volume!(ax2, field.x; algorithm = :mip, colormap = :viridis, transparency = true)
            end
        end
    end

    fig
end

function draw_cells_3D(sim; ax=nothing)
    scene = ax === nothing ? LScene(size=(800,800)) : ax

    model = sim.model
    state = model.state
    pxs = state.extra[:pixels_by_cell]
    colors = sim.cfg[:runtime][:cell_kind_color]

    for (cid, coords) in pxs
        clr = get(colors, sim.model.state.kinds[cid], RGBA(1, 0, 0, 1.0))
        for pos in coords
            x, y, z = Tuple(pos)
            mesh!(scene, Rect3f(Vec{3, Float32}(x, y, z), Vec{3, Float32}(1, 1, 1)), color=clr, transparency=true)
        end
    end

    scene
end

function prep_vis_sim(sim;
    burnin_steps = 20,
    show_3d = true,
    num_cell_kinds = nothing,
    show_activity = false,
    show_fields = false,
    cell_kind_color = default_colors,
    field_color = default_colors_field)

    burnin!(sim, burnin_steps)
    t = Observable(0.0)

    # Assign each cell ID a distinct color using DistinguishableColors

    cell_color_map = Dict{Int, RGBA{Float32}}()
    
    update_cell_colors!() = begin
        ids = [c.id for c in sim.cells]
        raw_colors = distinguishable_colors(length(ids))
        for (i, id) in enumerate(ids)
            rgb = raw_colors[i]
            cell_color_map[id] = RGBAf0(rgb.r, rgb.g, rgb.b, 0.9f0)
        end
    end
    

    update_cell_colors!()  # call once before the first step


    if show_3d
        fig = Figure(size = (1200, 800))  # Big canvas
        ax1 = Axis3(fig[1, 1], title = "Cells", aspect = :data)
        ax2 = Axis3(fig[1, 2], title = "Vaxes", aspect = :data)
        
        # Stat table goes in full row under the plots
        stats_ax = fig[2, 1:2] = GridLayout()
        
        # Create stat labels
        step_label = Label(stats_ax[1, 1], "Step: 0", halign = :left)
        count_label = Label(stats_ax[1, 2], "Cell count: 0", halign = :left)
        vol_label = Label(stats_ax[1, 3], "Avg vol: 0", halign = :left)
        max_vol_label = Label(stats_ax[1, 4], "Max vol: 0", halign = :left)
        receptor_label = Label(stats_ax[1, 5], "Avg receptor: 0", halign = :left)
        


        size_x, size_y, size_z = size(sim.model.grid)
        limits!(ax1, 0, size_x, 0, size_y, 0, size_z)

        # Track all mesh plots for replacement
        cell_meshes_obs = Observable(Plot[])

        # Field volumes
        volumes = Dict{Symbol, Any}()
        if show_fields
            for (key, v) in sim.cfg[:vaxes]
                if get(v, :show, false)
                    field = sim.vaxes[v[:id]]
                    vol = volume!(ax2, field.x; algorithm = :mip, colormap = :viridis, transparency = true)
                    volumes[key] = vol
                end
            end
        end

        on(t) do _
            time_macro_step!(sim)

            # Delete old meshes
            for p in cell_meshes_obs[]
                delete!(ax1.scene, p)
            end
            cell_meshes_obs[] = Plot[]

            update_cell_colors!()

            # Draw 3D cubes for cell voxels
            cell_voxels = get_stat(sim.model, :pixels_by_cell)
            for (cid, voxels) in cell_voxels
                kind = cell_kind(sim.model, cid)
                color = get(cell_color_map, cid, RGBAf0(1, 1, 1, 1))
                for v in voxels
                    x, y, z = Tuple(v)
                    p = mesh!(ax1, Rect3f(Vec3f(x, y, z), Vec3f(1, 1, 1)), color=color, transparency=true)
                    push!(cell_meshes_obs[], p)
                end
            end

            # Draw borders in black for each kind (from 2 onward)
            for k in 2:(num_cell_kinds === nothing ? maximum(sim.model.state.cell_to_kind) : num_cell_kinds)
                border_vox = draw_cell_borders(sim.model, k, RGBAf0(1, 1, 1, 1))  # returns RGBA array
                border_indices = findall(x -> x.alpha > 0.0, border_vox)
                for idx in border_indices
                    x, y, z = Tuple(idx)
                    b = mesh!(ax1, Rect3f(Vec3f(x, y, z), Vec3f(1, 1, 1)), color=RGBAf0(0, 0, 0, 1.0))
                    push!(cell_meshes_obs[], b)
                end
            end


            # Update fields
            if show_fields
                for (key, volplot) in volumes
                    field = sim.vaxes[sim.cfg[:vaxes][key][:id]]
                    vol_data = convert(Array{Float32, 3}, field.x)
                    volplot[4][] = vol_data  
                end
            end

            step_label.text[] = "Step: $(Int(round(t[])))"
            count_label.text[] = "Cell count: $(length(sim.cells))"
            
            pixels_by_cell = get_stat(sim.model, :pixels_by_cell)
            cell_volumes = [length(get(pixels_by_cell, c.id, [])) for c in sim.cells]
            
            vol_label.text[] = "Avg vol: $(round(mean(cell_volumes), digits=1))"
            max_vol_label.text[] = "Max vol: $(maximum(cell_volumes))"            
            
            receptor_vals = [sum(c.receptors) for c in sim.cells]
            receptor_label.text[] = "Avg receptor: $(round(mean(receptor_vals), digits=3))"
            
            


        end

    else
        fig = Figure()
        ax = fig[1, 1] = Axis(fig)
        vox = Observable(nothing)

        on(t) do _
            time_macro_step!(sim)
            vox[] = draw_cpm(sim,
                    num_cell_kinds = num_cell_kinds,
                    show_fields = show_fields,
                    show_activity = show_activity,
                    cell_kind_color = cell_kind_color,
                    field_color = field_color,
                    vox = vox[])
        end
    end

    return fig, t
end

function anim_sim(sim, t; num_of_steps = 1000, fps = 1.0/30., should_we_continue = always_true)
    for i in 1:num_of_steps
        t[] += 1
        if !should_we_continue(sim)
            break;
        end
        sleep(fps)
    end
end


function record_sim(sim, t, img; filename = "simulation.mp4", timestamps = 0.0:1.0/30.:33.0, should_we_continue = always_true)
    record(img, filename, timestamps; framerate = 1.0/timestamps.step.hi) do tt
        t[] += tt
        if !should_we_continue(sim)
            stop()
        end
    end
end

end # module

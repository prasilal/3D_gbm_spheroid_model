const state_counts_log = []

const step_counter = Ref(0)

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

function parse_node_label(label::String)
    m = match(r"(.*) at step (\d+)", label)  # <--- NOT parentheses
    if m === nothing
        error("Could not parse label: $label")
    end
    return m.captures[1], parse(Int, m.captures[2])
end



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

    # Now plot
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

            # 👇 Smart horizontal offset based on how many transitions happen at this step
            x_offset = (i - (n+1)/2) * 0.05

            text!(ax, mid_x - x_offset, mid_y + 0.02, text = string(value), align = (:center, :bottom), fontsize=15)
        end
    end

    display(fig)
end

function print_average_volumes_by_type(sim)
    volumes = Dict{Symbol, Vector{Int}}()

    for cell in values(sim.cells)
        ctype = cell.state[:type]
        vol = length(get(get_stat(sim.model, :pixels_by_cell), cell.id, []))
        print("$ctype: ")
        println(vol)
        push!(get!(volumes, ctype, []), vol)
    end

    println("📊 Average cell volumes by type:")
    for (ctype, vols) in volumes
        avg = mean(vols)
        println("  $ctype: $(round(avg, digits=2)) voxels")
    end
end

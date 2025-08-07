module CellLogic

function cur_cum_state(sim, the_cell, in_vaxes) # infers what the state should be 
    cfg = sim.cfg
    reactions = sim.irep[:reactions]
    id = the_cell.id

    cells_border_pixels = get_stat(sim.model, :border_pixels_by_cell)
    cell_border_list = get(get_stat(sim.model, :cell_neighbor_list), id, Dict())

    n = cfg[:vaxes] |> length
    out_vaxes = get_outer_vaxes(in_vaxes, get(cells_border_pixels, id, []))
    bound_vaxes = get_bound_vaxes(n, sim.cells, cell_border_list, cells_border_pixels)
    out_bound_vaxes = out_vaxes .+ bound_vaxes

    weights = zeros(Float64, n)
    for vax in 1:n
        rs = reactions_for(reactions, vax)
        wp = 0.0
        for r in rs
            if r[RW] == 0
                continue
            end

            lp = 1.0
            for (n, v) in r[RREACT]
                lp *= (v == vax ? the_cell.receptors[v] : out_bound_vaxes[v]) ^ n
            end
            wp += lp * r[RW]

            if r[RR_ABSORB]
                absorb_vaxes!(vax, out_vaxes, the_cell.receptors,
                              r, sim.vaxes, get(cells_border_pixels, id, []))
            end
        end
        weights[vax] = wp
    end

    idx = argmax(weights)

    if iszero(weights) || weights[idx] < cfg[:rule_graph][:min_weight] || weights[idx] == 0.0
        the_cell.state[:cum_state], the_cell.state[:cum_state_weight]
    else
        idx, weights[idx]
    end
end

function cell_divide!(cell, edge, sim) # cell division
    rs = cell.receptors |> sum
    max_rs = typeof(edge[3]) <: Number ? edge[3] : edge[3][1]
    max_rs_div = typeof(edge[3]) <: Number ? 0.5 : edge[3][2]

    if rs < max_rs
        return nothing
    end

    println("diving cell $(cell.id). num of receptors $rs threshold $(edge[3])")

    nid = divide_cell!(sim.model, cell.id)
    new_cell = Cell(nid, deepcopy(cell.receptors), deepcopy(cell.state))
    mult = fill(0.5, length(cell.receptors))

    for (r, m) in edge[4:end]
        mult[r] = m
    end

    cell.receptors .*= mult

    if sum(cell.receptors) >= max_rs
        cell.receptors .*= max_rs_div
    end

    @. mult = 1.0 - mult
    new_cell.receptors .*= mult

    if sum(new_cell.receptors) >= max_rs
        new_cell.receptors .*= max_rs_div
    end

    new_cell.state[:cum_state] = edge[1]

    new_cell
end

function necrosis!(cell, sim, vaxes) # spills vaxes from dead cells
    cellpixels = get(get_stat(sim.model, :pixels_by_cell), cell.id, [])

    for idx in 1:length(cell.receptors)
        produce_vax!(vaxes[idx], cellpixels, cell.receptors[idx])
    end
end

function hard_kill!(cell, sim) # deletes cell (apoptosis)
    kill_cell!(sim.model, cell.id)
end

function hard_kill_n_del!(cid, sim)
    cell = sim.cells[findfirst(c -> c.id == cid, sim.cells)]
    hard_kill!(cell, sim)
    deleteat!(sim.cells, findfirst(c -> c.id == cell.id, sim.cells))
end


function soft_kill!(cell, sim) # cell still occupies space, but is dead
    set_cell_kind!(sim.model, cell.id, CELL_TYPE_DEAD)
end

function cell_kill!(cell, edge, sim, vaxes) 
    if rand() > edge[3]
        return false
    end

    println("killing cell $(cell.id) by $edge")

    @match edge[4] begin
        :apoptosis => nothing
        :necrosis => necrosis!(cell, sim, vaxes)
    end

    @match edge[5] begin
        :hard => hard_kill!(cell, sim)
        :soft => soft_kill!(cell, sim)
    end

    true
end

function update_by_zg!(sim, cell, vaxes) # updates cell by zg
    killed, new_cell = false, nothing

    for edge in sim.irep[:zg][cell.state[:cum_state]]
        action = edge[2]
        if action == :prod_r
            produce_r_vax!(cell, edge)
        elseif action == :prod_v
            produce_v_vax!(cell, edge, sim, vaxes)
        elseif action == :divide
            new_cell = cell_divide!(cell, edge, sim)
        elseif action == :kill
            killed = cell_kill!(cell, edge, sim, vaxes)
        end
    end

    res = killed ? [] : [cell]
    if new_cell != nothing
        push!(res, new_cell)
    end

    res
end

function get_resting_time(graph, zg, cstate) # time until next state change
    zg_actions = zg[cstate]
    rt_pos = findfirst(x -> x[2] == :resting_time, zg_actions)

    if rt_pos == nothing
        return graph[:resting_time]
    end

    zg_actions[rt_pos][3]
end

function update_cum_state(sim, cell, in_vaxes) # updates cell state
    if cell.state[:resting_time] > 0
        cell.state[:resting_time] -= 1
        return
    end

    cstate, cweight = cur_cum_state(sim, cell, in_vaxes)
    println("Cstate = $(convert_idx_to_state(cstate, sim)) ")
    println("Cum_state = $(convert_idx_to_state(cell.state[:cum_state], sim)) ")
    if cell.state[:cum_state] != cstate
        cell.state[:resting_time] = get_resting_time(sim.cfg[:rule_graph], sim.irep[:zg], cstate)

        cell_state_update = get(sim.cfg, :cell_state_update, nothing)
        if cell_state_update != nothing
            for f in get(cell_state_update, :hooks, [])
                f(sim, cell, cstate, cweight, get(cell_state_update, :state, nothing))
            end
        end

        if get(get(sim.cfg[:runtime], :debug, Dict()), :show_cell_state_change, true) == true
            println("Update cum state for cell $(cell.id): old($(convert_idx_to_state(cell.state[:cum_state], sim)), $(cell.state[:cum_state_weight])); new($(convert_idx_to_state(cstate, sim)), $cweight)")
        end
    end
    cell.state[:cum_state] = cstate
    cell.state[:cum_state_weight] = cweight
    set_cell_kind!(sim.model, cell.id, cum_state_to_idx(cstate))
end

function update_r_vaxes(sim, cell) # updates receptors (decay)
    vaxes = sim.irep[:vaxes]

    for (idx, vinfo) in enumerate(vaxes)
        cell.receptors[idx] *= 1 - vinfo[:rd]
    end
end

end # module

function division_hook(cell, edge, sim)

    if cell.state[:phase] == :interphase
        if cell.state[:type] == :msc
            cell.state[:cell_cycle_timer] -= check_conditions(cell,sim, [:oxygen, :glucose, :lactate, :ecm_barrier, :crowding])
        else
            cell.state[:cell_cycle_timer] -= check_conditions(cell,sim)
        end
        
        if cell.state[:cell_cycle_timer] <=0
            if rand() < adjusted_division_probability(cell,sim)
                cell.state[:phase] = :mitosis
            else
                # quiescence?
                reset_timer!(cell, sim)
            end
        end
    end

    if cell.state[:phase] == :mitosis
        if sim.model.state.cell_volume[cell.id] <=2
            println("cell is too small to divide")
            cell.state[:phase] = :interphase
            reset_timer!(cell, sim)
            return nothing
        end
        nid = divide_cell!(sim.model, cell.id)
        new_cell = Cell(nid, deepcopy(cell.receptors), deepcopy(cell.state))
        new_cell.state[:cum_state] = edge[1] # u GSC možná jiné než mateřská buňka
        sim.cell_map[nid] = new_cell
        println("$(convert_idx_to_state(cell.state[:cum_state], sim)) divided")
        cell.state[:phase] = :interphase
        new_cell.state[:phase] = :interphase
        reset_timer!(cell, sim)
        reset_timer!(new_cell, sim)

        return new_cell
    end
    
end

function check_conditions(cell, sim, factors = [:oxygen, :glucose, :lactate, :ecm_barrier, :growth_signal, :crowding])

    env_score = get_environmental_score(cell, sim, factors)
    env_threshold = 0.8 # could be a parameter
    #println("Environmental score: $env_score")
    
    cfg = sim.cfg
	
    growth_signal = local_vax_concentration(sim, cell, :growth_signal)
    growth_speedup = sigmoid_boost(growth_signal, get(get(cfg[:vaxes], :growth_signal, Dict()), :thresholds, 5.0),3.0,0.0,1.0)

    ecm_barrier = local_vax_concentration(sim, cell, :ecm_barrier)
    barrier_slowdown = sigmoid_boost(ecm_barrier, get(get(cfg[:vaxes], :ecm_barrier, Dict()), :thresholds, 7.0), 3.0, 0.0,1.0)

    if env_score >= env_threshold
        return 1.0 + growth_speedup - barrier_slowdown
    elseif env_score <= 0.0
        return -1.0 + growth_speedup - barrier_slowdown
    else
        scaled = env_score / env_threshold
        slowdown = -1.0 + scaled  
        return Base.clamp(slowdown, -1.0, 1.0) + growth_speedup - barrier_slowdown
    end

end

function get_environmental_score(cell, sim, factors::Vector{Symbol})
    v = Dict{Symbol, Float64}()
    
    # Local concentrations
    local_conc = Dict(
        :oxygen => local_vax_concentration(sim, cell, :oxygen),
        :glucose => local_vax_concentration(sim, cell, :glucose),
        :lactate => local_vax_concentration(sim, cell, :lactate),
        :ecm_barrier => local_vax_concentration(sim, cell, :ecm_barrier),
        :growth_signal => local_vax_concentration(sim, cell, :growth_signal),
    )

    cfg = sim.cfg
    thresholds = Dict(
        :oxygen => get(get(cfg[:vaxes], :oxygen, Dict()), :thresholds, [0.05, 1.0]),
        :glucose => get(get(cfg[:vaxes], :glucose, Dict()), :thresholds, [0.074, 0.5]),
        :lactate => get(get(cfg[:vaxes], :lactate, Dict()), :thresholds, 3.0),
        :ecm_barrier => get(get(cfg[:vaxes], :ecm_barrier, Dict()), :thresholds, 6.0),
        :growth_signal => get(get(cfg[:vaxes], :growth_signal, Dict()), :thresholds, 4.0),
    )

    for f in factors
        if f == :oxygen
            lo, hi = thresholds[:oxygen]
            v[:oxygen] = Base.clamp((local_conc[:oxygen] - lo) / (hi - lo), 0.0, 1.0)
        elseif f == :glucose
            lo, hi = thresholds[:glucose]
            v[:glucose] = Base.clamp((local_conc[:glucose] - lo) / (hi - lo), 0.0, 1.0)
        elseif f == :lactate
            max_l = thresholds[:lactate]
            v[:lactate] = Base.clamp(1 - (local_conc[:lactate] / max_l)^2, 0.0, 1.0)
        elseif f == :ecm_barrier
            max_e = thresholds[:ecm_barrier]
            v[:ecm_barrier] = Base.clamp(1 - (local_conc[:ecm_barrier] / max_e)^1.5, 0.0, 1.0)
        elseif f == :growth_signal
            max_g = thresholds[:growth_signal]
            v[:growth_signal] = sigmoid_boost(local_conc[:growth_signal], max_g)
        elseif f == :crowding
            # Get neighbor info
            cell_border_list = get(get_stat(sim.model, :cell_neighbor_list), cell.id, Dict())
            neighbor_count = count(id -> id > CELL_TYPE_DEAD, keys(cell_border_list))
            total_contact_pixels = sum((v for (id, v) in cell_border_list if id > CELL_TYPE_DEAD); init=0.0)
            total_border_pixels = sum((values(cell_border_list)); init=0.0)
            contact_fraction = total_contact_pixels / max(total_border_pixels, 1)
            # Neighbor count clamp (6 → 1.0, 8+ → 0.0)
            penalty_neighbors = Base.clamp((8 - neighbor_count) / (8 - 6), 0.0, 1.0)
            # Contact coverage clamp (≤0.8 → 1.0, 1.0 → 0.0)
            penalty_contact = Base.clamp((1.0 - contact_fraction) / 0.2, 0.0, 1.0)

            v[:crowding] = 0.5 * penalty_neighbors + 0.5 * penalty_contact

        end
    end

    used = [v[f] for f in factors if haskey(v, f)]

    if :growth_signal in factors
        gboost = v[:growth_signal]
        base = [v[f] for f in factors if f != :growth_signal]
        score = Base.clamp(prod(base)^(1 / length(base)) * gboost, 0.0, 1.0)
    else
        score = Base.clamp(prod(used)^(1 / length(used)), 0.0, 1.0)
    end

    return score
end

function sigmoid_boost(g, midpoint, steepness=4.0, base=1.0, max_boost=0.5)
    return base + max_boost / (1 + exp(-steepness * (g - midpoint)))
end


function adjusted_division_probability(cell, sim)

    cfg = sim.cfg
    cum_state = convert_idx_to_state(cell.state[:cum_state], sim)
    baseline_probability = cfg[:rule_graph][:zg][cum_state][findfirst(e -> e[2] == :divide, cfg[:rule_graph][:zg][cum_state])][3]

    env_score = get_environmental_score(cell, sim, [:oxygen, :glucose, :lactate, :ecm_barrier, :growth_signal, :crowding])

    shift = 0.1 * (env_score - 0.5)  # ∈ [-0.05, +0.05]
    adjusted_probability = Base.clamp(baseline_probability + shift, 0.0, 1.0)
    return adjusted_probability

end

function reset_timer!(cell, sim)
    cfg = sim.cfg
    interphase_duration = get(get(cfg[:vaxes], convert_idx_to_state(cell.state[:cum_state], sim), Dict()), :interphase_duration, 2)  # TODO: change the default parameter to something biologically meaningful
    cell.state[:cell_cycle_timer] = interphase_duration
end

"""
division_hook(cell, edge, sim)

Hook function that governs the logic for cell division. Evaluates whether a cell in interphase 
progresses toward mitosis based on environmental conditions, and if it divides, spawns a new cell 
with copied state and adjusted `cum_state`.

Arguments:
- `cell`: the current cell object.
- `edge`: tuple containing information about which state to assign to the daughter cell.
- `sim`: the simulation object.

Returns:
- A new `Cell` object if division occurs, otherwise `nothing`.
"""
function division_hook(cell, edge, sim)

    # === INTERPHASE PHASE ===
    if cell.state[:phase] == :interphase
        # MSC division depends on additional local signals
        if cell.state[:type] == :msc
            cell.state[:cell_cycle_timer] -= check_conditions(cell, sim, [:oxygen, :glucose, :lactate, :ecm_barrier, :crowding])
        else
            cell.state[:cell_cycle_timer] -= check_conditions(cell, sim)
        end

        # Cell attempts to divide once timer expires
        if cell.state[:cell_cycle_timer] <= 0
            if rand() < adjusted_division_probability(cell, sim)
                cell.state[:phase] = :mitosis
            else
                reset_timer!(cell, sim)  # Reset if division fails
            end
        end
    end

    # === MITOSIS PHASE ===
    if cell.state[:phase] == :mitosis
        if sim.model.state.cell_volume[cell.id] <= 2
            println("cell is too small to divide")
            cell.state[:phase] = :interphase
            reset_timer!(cell, sim)
            return nothing
        end

        # Perform division
        nid = divide_cell!(sim.model, cell.id)
        new_cell = Cell(nid, deepcopy(cell.receptors), deepcopy(cell.state))
        new_cell.state[:cum_state] = edge[1]  # Set daughter state from edge
        sim.cell_map[nid] = new_cell

        println("$(convert_idx_to_state(cell.state[:cum_state], sim)) divided")

        # Reset phase/timer for both cells
        cell.state[:phase] = :interphase
        new_cell.state[:phase] = :interphase
        reset_timer!(cell, sim)
        reset_timer!(new_cell, sim)

        return new_cell
    end
end


"""
check_conditions(cell, sim, factors)

Returns an adjustment to the cell cycle timer based on local environmental conditions.
Positive return values speed up progression toward mitosis; negative values delay it.

- `factors`: environmental vaxes to consider (oxygen, glucose, lactate, etc.)
"""
function check_conditions(cell, sim, factors = [:oxygen, :glucose, :lactate, :ecm_barrier, :growth_signal, :crowding])
    env_score = get_environmental_score(cell, sim, factors)
    env_threshold = 0.8  # Above this, environment is considered favorable

    # Growth-promoting factors
    growth_signal = local_vax_concentration(sim, cell, :growth_signal)
    growth_speedup = sigmoid_boost(growth_signal, get(get(sim.cfg[:vaxes], :growth_signal, Dict()), :thresholds, 5.0), 3.0, 0.0, 1.0)

    # ECM barriers slow growth
    ecm_barrier = local_vax_concentration(sim, cell, :ecm_barrier)
    barrier_slowdown = sigmoid_boost(ecm_barrier, get(get(sim.cfg[:vaxes], :ecm_barrier, Dict()), :thresholds, 7.0), 3.0, 0.0, 1.0)

    # Timer effect calculation
    if env_score >= env_threshold
        return 1.0 + growth_speedup - barrier_slowdown
    elseif env_score <= 0.0
        return -1.0 + growth_speedup - barrier_slowdown
    else
        scaled = env_score / env_threshold
        slowdown = -1.0 + scaled
        return clamp(slowdown, -1.0, 1.0) + growth_speedup - barrier_slowdown
    end
end


"""
get_environmental_score(cell, sim, factors)

Combines normalized local vax concentrations into a single score [0.0, 1.0].
High score indicates good growth conditions.

- `factors`: list of vaxes to include in scoring (e.g., :oxygen, :crowding)
"""
function get_environmental_score(cell, sim, factors::Vector{Symbol})
    v = Dict{Symbol, Float64}()

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
            v[:oxygen] = clamp((local_conc[:oxygen] - lo) / (hi - lo), 0.0, 1.0)
        elseif f == :glucose
            lo, hi = thresholds[:glucose]
            v[:glucose] = clamp((local_conc[:glucose] - lo) / (hi - lo), 0.0, 1.0)
        elseif f == :lactate
            max_l = thresholds[:lactate]
            v[:lactate] = clamp(1 - (local_conc[:lactate] / max_l)^2, 0.0, 1.0)
        elseif f == :ecm_barrier
            max_e = thresholds[:ecm_barrier]
            v[:ecm_barrier] = clamp(1 - (local_conc[:ecm_barrier] / max_e)^1.5, 0.0, 1.0)
        elseif f == :growth_signal
            max_g = thresholds[:growth_signal]
            v[:growth_signal] = sigmoid_boost(local_conc[:growth_signal], max_g)
        elseif f == :crowding
            # Crowding penalizes dense environments
            nb_list = get(get_stat(sim.model, :cell_neighbor_list), cell.id, Dict())
            neighbor_count = count(id -> id > CELL_TYPE_DEAD, keys(nb_list))
            total_contact = sum(v for (id, v) in nb_list if id > CELL_TYPE_DEAD)
            total_border = sum(values(nb_list))
            contact_fraction = total_contact / max(total_border, 1)

            penalty_neighbors = clamp((8 - neighbor_count) / 2, 0.0, 1.0)
            penalty_contact = clamp((1.0 - contact_fraction) / 0.2, 0.0, 1.0)

            v[:crowding] = 0.5 * penalty_neighbors + 0.5 * penalty_contact
        end
    end

    used = [v[f] for f in factors if haskey(v, f)]

    # Aggregate into a single score using geometric mean + growth boost
    if :growth_signal in factors
        base = [v[f] for f in factors if f != :growth_signal]
        return clamp(prod(base)^(1 / length(base)) * v[:growth_signal], 0.0, 1.0)
    else
        return clamp(prod(used)^(1 / length(used)), 0.0, 1.0)
    end
end


"""
sigmoid_boost(g, midpoint, steepness=4.0, base=1.0, max_boost=0.5)

Returns a sigmoid-modulated boost based on input `g`, used to scale division and timer effects.
"""
function sigmoid_boost(g, midpoint, steepness=4.0, base=1.0, max_boost=0.5)
    return base + max_boost / (1 + exp(-steepness * (g - midpoint)))
end

"""
adjusted_division_probability(cell, sim)

Computes final probability that a cell divides after interphase, based on environment score.

Returns:
- Float in [0.0, 1.0]
"""
function adjusted_division_probability(cell, sim)
    cfg = sim.cfg
    cum_state = convert_idx_to_state(cell.state[:cum_state], sim)
    base_prob = cfg[:rule_graph][:zg][cum_state][findfirst(e -> e[2] == :divide, cfg[:rule_graph][:zg][cum_state])][3]

    env_score = get_environmental_score(cell, sim, [:oxygen, :glucose, :lactate, :ecm_barrier, :growth_signal, :crowding])
    shift = 0.1 * (env_score - 0.5)  # Scale between -0.05 to +0.05
    return clamp(base_prob + shift, 0.0, 1.0)
end

"""
reset_timer!(cell, sim)

Resets the cell’s cell cycle timer to its configured `interphase_duration`.
"""
function reset_timer!(cell, sim)
    cfg = sim.cfg
    cum_state = convert_idx_to_state(cell.state[:cum_state], sim)
    interphase_duration = get(get(cfg[:vaxes], cum_state, Dict()), :interphase_duration, 2)
    cell.state[:cell_cycle_timer] = interphase_duration
end


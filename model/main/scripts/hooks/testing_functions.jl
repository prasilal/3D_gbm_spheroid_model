const TRANSITION_FN_MAP = Dict{Symbol, Function}()

const transition_log = []

# GBM
TRANSITION_FN_MAP[:gbm] = function (sim, cell, current_state_idx, weight, _)
    vaxmap = sim.irep[:vaxmap]
    current_state = convert_idx_to_state(current_state_idx, sim)

    # Receptor values
    o = local_vax_concentration(sim, cell, :oxygen)
    g = local_vax_concentration(sim, cell, :glucose)
    if current_state !== nothing
        return
    end

    if current_state in [:gbm_p, :gbm_q]
        if o < 0.0001 || g < 0.0001
            cell.state[:cum_state] = vaxmap[:gbm_n]
        elseif o < 0.0002 || g < 0.0002
            cell.state[:cum_state] = vaxmap[:gbm_q]
        else
            cell.state[:cum_state] = vaxmap[:gbm_p]

        end
        log_transition(cell.id, current_state_idx, cell.state[:cum_state], sim)
    end

    cell.receptors .*= 0.5
    cell.receptors[vaxmap[convert_idx_to_state(cell.state[:cum_state], sim)]] = 1.0

    cell.state[:cum_state_weight] = weight
    set_cell_kind!(sim.model, cell.id, cum_state_to_idx(cell.state[:cum_state]))
end

# MSC
TRANSITION_FN_MAP[:msc] = function (sim, cell, current_state_idx, weight, _)
    vaxmap = sim.irep[:vaxmap]
    current_state = convert_idx_to_state(current_state_idx, sim)

    # Receptor values
    o = local_vax_concentration(sim, cell, :oxygen)
    g = local_vax_concentration(sim, cell, :glucose)
    if current_state !== nothing
        return
    end

    if current_state in [:gbm_p, :gbm_q]
        if o < 0.0001 || g < 0.0001
            cell.state[:cum_state] = vaxmap[:gbm_n]
        elseif o < 0.0002 || g < 0.0002
            cell.state[:cum_state] = vaxmap[:gbm_q]
        else
            cell.state[:cum_state] = vaxmap[:gbm_p]

        end
        log_transition(cell.id, current_state_idx, cell.state[:cum_state], sim)
    end

    cell.receptors .*= 0.5
    cell.receptors[vaxmap[convert_idx_to_state(cell.state[:cum_state], sim)]] = 1.0

    cell.state[:cum_state_weight] = weight
    set_cell_kind!(sim.model, cell.id, cum_state_to_idx(cell.state[:cum_state]))
end

# GSC
TRANSITION_FN_MAP[:gsc] = function (sim, cell, current_state_idx, weight, _)
    vaxmap = sim.irep[:vaxmap]
    current_state = convert_idx_to_state(current_state_idx, sim)

    # Receptor values
    o = local_vax_concentration(sim, cell, :oxygen)
    g = local_vax_concentration(sim, cell, :glucose)
    if current_state !== nothing
        return
    end

    if current_state in [:gbm_p, :gbm_q]
        if o < 0.0001 || g < 0.0001
            cell.state[:cum_state] = vaxmap[:gbm_n]
        elseif o < 0.0002 || g < 0.0002
            cell.state[:cum_state] = vaxmap[:gbm_q]
        else
            cell.state[:cum_state] = vaxmap[:gbm_p]

        end
        log_transition(cell.id, current_state_idx, cell.state[:cum_state], sim)
    end

    cell.receptors .*= 0.5
    cell.receptors[vaxmap[convert_idx_to_state(cell.state[:cum_state], sim)]] = 1.0

    cell.state[:cum_state_weight] = weight
    set_cell_kind!(sim.model, cell.id, cum_state_to_idx(cell.state[:cum_state]))
end

function test_transition_dispatcher(sim, cell, in_vaxes)
    cell_type = cell.state[:type]

    transition_fn = get(TRANSITION_FN_MAP, cell_type, nothing)
    if transition_fn === nothing
        println("No transition function for cell type: $cell_type")
        return
    end

    current_state_idx = cell.state[:cum_state]
    absorb_vaxes_for_cell!(sim, cell, in_vaxes)

    transition_fn(sim, cell, current_state_idx, 1.0, nothing)
    log_cell_states(sim, cell, current_state_idx, 1.0, nothing)
end

function test_division_hook(cell, edge, sim)

    if cell.state[:phase] == :interphase
        cell.state[:cell_cycle_timer] -= test_check_conditions(cell,sim)

        if cell.state[:cell_cycle_timer] <=0
            if rand() < test_adjusted_division_probability(cell,sim)
                cell.state[:phase] = :mitosis
            else
                test_reset_timer!(cell, sim)
            end
        end
    end

    if cell.state[:phase] == :mitosis
        if sim.model.state.cell_volume[cell.id] <=2
            println("cell is too small to divide")
            cell.state[:phase] = :interphase
            test_reset_timer!(cell, sim)
            return nothing
        end
        nid = divide_cell!(sim.model, cell.id)
        new_cell = Cell(nid, deepcopy(cell.receptors), deepcopy(cell.state))
        new_cell.state[:cum_state] = edge[1] 
        println("$(convert_idx_to_state(cell.state[:cum_state], sim)) divided")
        cell.state[:phase] = :interphase
        new_cell.state[:phase] = :interphase
        test_reset_timer!(cell, sim)
        test_reset_timer!(new_cell, sim)

        return new_cell
    end
    
end

function test_check_conditions(cell, sim)

    env_score = test_get_environmental_score(cell, sim)
    env_threshold = 0.8 # could be a parameter

    if env_score >= env_threshold
        return 1.0
    elseif env_score <= 0.0
        return -1.0
    else
        scaled = env_score / env_threshold
        slowdown = -1.0 + scaled  
        return Base.clamp(slowdown, -1.0, 0.0)
    end

end

function test_get_environmental_score(cell, sim)
    o = local_vax_concentration(sim, cell, :oxygen)
    g = local_vax_concentration(sim, cell, :glucose)

    # Nutrient conditions
    cfg = sim.cfg
    o_thresholds = get(get(cfg[:vaxes], :oxygen, Dict()), :thresholds, [0.01, 0.065])
    g_thresholds = get(get(cfg[:vaxes], :glucose, Dict()), :thresholds, [0.5, 5.0])

    o_norm = Base.clamp((o - o_thresholds[1]) / (o_thresholds[2] - o_thresholds[1]), 0.0, 1.0)
    g_norm = Base.clamp((g - g_thresholds[1]) / (g_thresholds[2] - g_thresholds[1]), 0.0, 1.0)

    # Mechanical stress
    cum_state = convert_idx_to_state(cell.state[:cum_state], sim)
    current_vol = sim.model.state.cell_volume[cell.id]
    target_vol = cfg[:rule_graph][:zg][cum_state][findfirst(e -> e[2] == :volume, cfg[:rule_graph][:zg][cum_state])][4]

    vol_ratio = current_vol / target_vol
    v_norm = Base.clamp((vol_ratio - 0.5) / 0.5, 0.0, 1.0) # TODO: opravit logicky vol_ratio
    

    env_score = (o_norm * g_norm * v_norm)^(1/3)
    return env_score
end

function test_adjusted_division_probability(cell, sim)

    cfg = sim.cfg
    cum_state = convert_idx_to_state(cell.state[:cum_state], sim)
    baseline_probability = cfg[:rule_graph][:zg][cum_state][findfirst(e -> e[2] == :divide, cfg[:rule_graph][:zg][cum_state])][3]

    env_score = test_get_environmental_score(cell, sim)

    shift = 0.1 * (env_score - 0.5)  # ∈ [-0.05, +0.05]
    adjusted_probability = Base.clamp(baseline_probability + shift, 0.0, 1.0)
    return adjusted_probability

end

function test_reset_timer!(cell, sim)
    cfg = sim.cfg
    interphase_duration = get(get(cfg[:vaxes], convert_idx_to_state(cell.state[:cum_state], sim), Dict()), :interphase_duration, 20)  # TODO: change the default parameter to something biologically meaningful
    cell.state[:cell_cycle_timer] = interphase_duration
end
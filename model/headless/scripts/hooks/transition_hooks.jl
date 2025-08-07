const TRANSITION_FN_MAP = Dict{Symbol, Function}()

const transition_log = []


function try_transition(sim, cell, prob, target_state, vaxmap; stress_bump = 0)
    if rand() < prob
        # log_transition(cell.id, cell.state[:cum_state], vaxmap[target_state], sim)
        cell.state[:cum_state] = vaxmap[target_state]
        cell.state[:resting_tracker] = sim.cfg[:rule_graph][:resting_time]
        if stress_bump > 0
            cell.state[:stress] = max(cell.state[:stress] - stress_bump, 0)
        end
    end
end


# GBM
TRANSITION_FN_MAP[:gbm] = function (sim, cell, current_state_idx, weight, _)
    vaxmap = sim.irep[:vaxmap]
    current_state = convert_idx_to_state(current_state_idx, sim)
    thresholds = sim.cfg[:cell_state_update][:thresholds]
    transitions = sim.cfg[:cell_state_update][:transitions]


    # Local vax concentration values
    o = local_vax_concentration(sim, cell, :oxygen)
    g = local_vax_concentration(sim, cell, :glucose)
    lactate = local_vax_concentration(sim, cell, :lactate)
    ecm = local_vax_concentration(sim, cell, :ecm_barrier)
    growth_signal = local_vax_concentration(sim, cell, :growth_signal)
    ecm_gsc = local_vax_concentration(sim, cell, :ecm_gsc)



    hypoxia = o < thresholds[:hypoxia]
    severe_hypoxia = o < thresholds[:severe_hypoxia]
    starvation = g < thresholds[:starvation]
    severe_starvation = g < thresholds[:severe_starvation]

    growth_threshold = thresholds[:growth_signal]
    g_boost = sigmoid_boost(growth_signal, growth_threshold, 3.0, 0.0, 0.2)

    stemness_threshold = thresholds[:ecm_gsc]
    stemness = sigmoid_boost(ecm_gsc, stemness_threshold, 3.0, 0.0, 0.8)


    env_score = get_environmental_score(cell, sim, [:oxygen, :glucose, :lactate, :ecm_barrier])
    if env_score < thresholds[:env_score_threshold][:high] || hypoxia || starvation
        cell.state[:stress] += 1
    else
        cell.state[:stress] = max(0, cell.state[:stress] - 1)
    end

    # --------------------
    # GBM Proliferative
    # --------------------

    if current_state == :gbm_p
        if severe_hypoxia
            if severe_starvation || cell.state[:stress] >= transitions[:gbm_stress_n] 
                try_transition(sim, cell, transitions[:gbm_p_to_n_prob], :gbm_n, vaxmap)
            else
                try_transition(sim, cell, transitions[:gbm_p_to_q_prob][:high]-0.05-g_boost, :gbm_q, vaxmap)
            end
        elseif hypoxia 
            if severe_starvation || cell.state[:stress] >= transitions[:gbm_stress_n] 
                try_transition(sim, cell, transitions[:gbm_p_to_q_prob][:high]-g_boost, :gbm_q, vaxmap)
            elseif starvation || cell.state[:stress] >= transitions[:gbm_stress_q] && env_score < thresholds[:env_score_threshold][:low] 
                try_transition(sim, cell, transitions[:gbm_p_to_q_prob][:medium]-g_boost, :gbm_q, vaxmap)

            else
                if rand() < stemness
                    try_transition(sim, cell, transitions[:gbm_to_gsc_prob], :gsc_p, vaxmap)
                    println("Transitioning to GSC Proliferative due to stemness: $stemness")
                elseif cell.state[:stress] >= transitions[:gbm_stress_q]
                    try_transition(sim, cell, transitions[:gbm_p_to_q_prob][:low]-g_boost, :gbm_q, vaxmap)
                end
                
            end
        elseif starvation
            if severe_starvation && cell.state[:stress] >= transitions[:gbm_stress_n] 
                try_transition(sim, cell, transitions[:gbm_p_to_q_prob][:medium]-g_boost, :gbm_q, vaxmap)
            elseif cell.state[:stress] >= transitions[:gbm_stress_q] && (env_score < thresholds[:env_score_threshold][:low] || severe_starvation) 
                try_transition(sim, cell, transitions[:gbm_p_to_q_prob][:low]-g_boost, :gbm_q, vaxmap)
            else
                if cell.state[:stress] >= transitions[:gbm_stress_q]
                    try_transition(sim, cell, transitions[:gbm_p_to_q_prob][:verylow]-g_boost, :gbm_q, vaxmap)
                end
                
            end
        end

    end

    # --------------------
    # GBM Quiescent
    # --------------------

    if current_state == :gbm_q
        if severe_hypoxia
            if severe_starvation || cell.state[:stress] >= transitions[:gbm_stress_n] 
                try_transition(sim, cell, transitions[:gbm_q_to_n_prob][:high], :gbm_n, vaxmap)
            elseif cell.state[:stress] >= transitions[:gbm_stress_q]
                if starvation || env_score < thresholds[:env_score_threshold][:low]
                    try_transition(sim, cell, transitions[:gbm_q_to_n_prob][:medium], :gbm_n, vaxmap)
                else
                    try_transition(sim, cell, transitions[:gbm_q_to_n_prob][:low], :gbm_n, vaxmap)
                end
            end
        elseif severe_starvation
            if hypoxia && cell.state[:stress] >= transitions[:gbm_stress_n]
                try_transition(sim, cell, transitions[:gbm_q_to_n_prob][:high], :gbm_n, vaxmap)
            elseif cell.state[:stress] >= transitions[:gbm_stress_n] && env_score < thresholds[:env_score_threshold][:low]
                try_transition(sim, cell, transitions[:gbm_q_to_n_prob][:low], :gbm_n, vaxmap)
            end   
        elseif !hypoxia
            if starvation || cell.state[:stress] >= transitions[:gbm_stress_q]
                try_transition(sim, cell, transitions[:gbm_q_to_p_prob][:low]+g_boost, :gbm_p, vaxmap, stress_bump=transitions[:gbm_q_to_p_stress_bump])
            elseif cell.state[:stress] >= transitions[:gbm_stress_q]÷2 || env_score < thresholds[:env_score_threshold][:low]
                try_transition(sim, cell, transitions[:gbm_q_to_p_prob][:medium]+g_boost, :gbm_p, vaxmap, stress_bump=transitions[:gbm_q_to_p_stress_bump])
            else 
                try_transition(sim, cell, transitions[:gbm_q_to_p_prob][:high]+g_boost, :gbm_p, vaxmap, stress_bump=transitions[:gbm_q_to_p_stress_bump])
            end
        end
    end

    cell.receptors .*= 0.5
    cell.receptors[vaxmap[convert_idx_to_state(cell.state[:cum_state], sim)]] = 1.0

    cell.state[:cum_state_weight] = weight
    set_cell_kind!(sim.model, cell.id, cum_state_to_idx(cell.state[:cum_state]))

end

function contact_fraction_prob(sim, cell)
    cell_border_list = get(get_stat(sim.model, :cell_neighbor_list), cell.id, Dict())
    total_border_pixels = sum(values(cell_border_list); init=0.0)

    contact_pixels = sum(
        (v for (id, v) in cell_border_list if id > CELL_TYPE_DEAD &&
         convert_idx_to_state(sim.cell_map[id].state[:cum_state], sim) in (:gbm_p, :gsc_p)); 
        init=0.0
    )

    return contact_pixels / max(total_border_pixels,1.0)
end


# MSC
TRANSITION_FN_MAP[:msc] = function (sim, cell, current_state_idx, weight, _)
    vaxmap = sim.irep[:vaxmap]
    current_state = convert_idx_to_state(current_state_idx, sim)
    thresholds = sim.cfg[:cell_state_update][:thresholds]
    transitions = sim.cfg[:cell_state_update][:transitions]

    # Local vax concentration values
    o = local_vax_concentration(sim, cell, :oxygen)
    lactate = local_vax_concentration(sim, cell, :lactate)

    hypoxia = o < thresholds[:hypoxia]
    severe_hypoxia = o < thresholds[:severe_hypoxia]
    acidity = sigmoid_boost(lactate, thresholds[:acidity], 3.0, 0.0, 1.0)

    gbm_contact_prob = contact_fraction_prob(sim, cell)


    env_score = get_environmental_score(cell, sim, [:oxygen, :glucose, :lactate, :ecm_barrier])
    if env_score < thresholds[:env_score_threshold][:medium] || hypoxia || gbm_contact_prob > transitions[:gbm_contact_prob][:medium]
        cell.state[:stress] += 1
    else
        cell.state[:stress] = max(0, cell.state[:stress] - 1)
    end

    # --------------------
    # MSC Pro-tumor
    # --------------------

    if current_state == :msc_p
        if !hypoxia
            if acidity < transitions[:acidity_thresholds][:low] || gbm_contact_prob < transitions[:gbm_contact_prob][:verylow]
                try_transition(sim, cell, transitions[:msc_p_to_anti_prob][:high], :msc_anti, vaxmap, stress_bump=transitions[:msc_p_to_anti_stress_bump])
            elseif acidity < transitions[:acidity_thresholds][:medium] || gbm_contact_prob < transitions[:gbm_contact_prob][:low] 
                try_transition(sim, cell, transitions[:msc_p_to_anti_prob][:medium], :msc_anti, vaxmap, stress_bump=transitions[:msc_p_to_anti_stress_bump])
            else
                try_transition(sim, cell, transitions[:msc_p_to_anti_prob][:verylow], :msc_anti, vaxmap, stress_bump=transitions[:msc_p_to_anti_stress_bump])
            end
        elseif acidity < transitions[:acidity_thresholds][:low] || gbm_contact_prob < transitions[:gbm_contact_prob][:verylow]
            try_transition(sim, cell, transitions[:msc_p_to_anti_prob][:low], :msc_anti, vaxmap, stress_bump=transitions[:msc_p_to_anti_stress_bump])
        end

    end

    # --------------------
    # MSC Anti-tumor
    # --------------------

    if current_state == :msc_anti
        if severe_hypoxia
            if cell.state[:stress] >= transitions[:msc_stress_n] || gbm_contact_prob > transitions[:gbm_contact_prob][:high] || acidity > transitions[:acidity_thresholds][:high] 
                try_transition(sim, cell, transitions[:msc_anti_to_p_prob][:high], :msc_p, vaxmap)
            elseif cell.state[:stress] >= transitions[:msc_stress_q] || gbm_contact_prob > transitions[:gbm_contact_prob][:medium] 
                try_transition(sim, cell, transitions[:msc_anti_to_p_prob][:medium], :msc_p, vaxmap)
            else
                try_transition(sim, cell, transitions[:msc_anti_to_p_prob][:low], :msc_p, vaxmap)
            end
        elseif acidity > transitions[:acidity_thresholds][:medium] || gbm_contact_prob > transitions[:gbm_contact_prob][:medium] 
            try_transition(sim, cell, transitions[:msc_anti_to_p_prob][:low], :msc_p, vaxmap)
        end
    end

    # --------------------
    # MSC Before spheroid formation
    # --------------------

    if current_state == :msc_before
        try_transition(sim, cell, 0.8, :msc_p, vaxmap)
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
    thresholds = sim.cfg[:cell_state_update][:thresholds]
    transitions = sim.cfg[:cell_state_update][:transitions]

    # Local vax concentration values
    o = local_vax_concentration(sim, cell, :oxygen)
    g = local_vax_concentration(sim, cell, :glucose)
    lactate = local_vax_concentration(sim, cell, :lactate)
    ecm = local_vax_concentration(sim, cell, :ecm_barrier)
    growth_signal = local_vax_concentration(sim, cell, :growth_signal)
    ecm_gsc = local_vax_concentration(sim, cell, :ecm_gsc)



    hypoxia = o < thresholds[:hypoxia]
    severe_hypoxia = o < thresholds[:severe_hypoxia]
    starvation = g < thresholds[:starvation]
    severe_starvation = g < thresholds[:severe_starvation]

    growth_threshold = thresholds[:growth_signal]
    g_boost = sigmoid_boost(growth_signal, growth_threshold, 3.0, 0.0, 0.2)

    stemness_threshold = thresholds[:ecm_gsc]
    stemness = sigmoid_boost(ecm_gsc, stemness_threshold, 3.0, 0.0, 0.2)


    env_score = get_environmental_score(cell, sim, [:oxygen, :glucose, :lactate, :ecm_barrier])
    if env_score < thresholds[:env_score_threshold][:high] || hypoxia || starvation
        cell.state[:stress] += 1
    else
        cell.state[:stress] = max(0, cell.state[:stress] - 1)
    end

    # --------------------
    # GSC Proliferative
    # --------------------

    if current_state == :gsc_p
        if severe_hypoxia
            if severe_starvation || cell.state[:stress] >= transitions[:gsc_stress_q] 
                try_transition(sim, cell, transitions[:gsc_p_to_gbm_prob]-stemness, :gbm_q, vaxmap)
            else
                try_transition(sim, cell, transitions[:gsc_p_to_q_prob][:higher]-g_boost, :gsc_q, vaxmap)
            end
        elseif hypoxia 
            if severe_starvation || cell.state[:stress] >= transitions[:gsc_stress_n] 
                try_transition(sim, cell, transitions[:gsc_p_to_q_prob][:higher]-g_boost, :gsc_q, vaxmap)
            elseif starvation || cell.state[:stress] >= transitions[:gsc_stress_q]  && env_score < thresholds[:env_score_threshold][:low] 
                try_transition(sim, cell, transitions[:gsc_p_to_q_prob][:high]-g_boost, :gsc_q, vaxmap)

            else
                if cell.state[:stress] >= transitions[:gsc_stress_q]
                    try_transition(sim, cell, transitions[:gsc_p_to_q_prob][:low]-g_boost, :gsc_q, vaxmap)
                end
                
            end
        elseif starvation
            if severe_starvation && cell.state[:stress] >= transitions[:gsc_stress_n] 
                try_transition(sim, cell, transitions[:gsc_p_to_q_prob][:veryhigh]-g_boost, :gsc_q, vaxmap)
            elseif cell.state[:stress] >= transitions[:gsc_stress_q] && (env_score < thresholds[:env_score_threshold][:low] || severe_starvation) 
                try_transition(sim, cell, transitions[:gsc_p_to_q_prob][:medium]-g_boost, :gsc_q, vaxmap)
            else
                if cell.state[:stress] >= transitions[:gsc_stress_n]
                    try_transition(sim, cell, transitions[:gsc_p_to_q_prob][:verylow]-g_boost, :gsc_q, vaxmap)
                end
                
            end
        end

    end

    # --------------------
    # GSC Quiescent
    # --------------------

    if current_state == :gsc_q
        if severe_hypoxia
            if severe_starvation || cell.state[:stress] >= transitions[:gsc_stress_n] 
                try_transition(sim, cell, transitions[:gsc_q_to_gbm_prob][:high]-stemness, :gbm_q, vaxmap)
            elseif cell.state[:stress] >= transitions[:gsc_stress_q]
                if starvation || env_score < thresholds[:env_score_threshold][:low]
                    try_transition(sim, cell, transitions[:gsc_q_to_gbm_prob][:medium]-stemness, :gbm_q, vaxmap)
                else
                    try_transition(sim, cell, transitions[:gsc_q_to_gbm_prob][:low]-stemness, :gbm_q, vaxmap)
                end
            end
        elseif severe_starvation
            if hypoxia && cell.state[:stress] >= transitions[:gsc_stress_n] 
                try_transition(sim, cell, transitions[:gsc_q_to_gbm_prob][:high]-stemness, :gbm_q, vaxmap)
            elseif cell.state[:stress] >= transitions[:gsc_stress_n]  && env_score < thresholds[:env_score_threshold][:low]
                try_transition(sim, cell, transitions[:gsc_q_to_gbm_prob][:low]-stemness, :gbm_q, vaxmap)
            end   
        elseif !hypoxia
            if starvation || cell.state[:stress] >= transitions[:gsc_stress_q] 
                try_transition(sim, cell, transitions[:gsc_q_to_p_prob][:low]+g_boost, :gsc_p, vaxmap, stress_bump=transitions[:gsc_q_to_p_stress_bump])
            elseif cell.state[:stress] >= transitions[:gsc_stress_q]  || env_score < thresholds[:env_score_threshold][:low] 
                try_transition(sim, cell, transitions[:gsc_q_to_p_prob][:medium]+g_boost, :gsc_p, vaxmap, stress_bump=transitions[:gsc_q_to_p_stress_bump])
            else 
                try_transition(sim, cell, transitions[:gsc_q_to_p_prob][:high]+g_boost, :gsc_p, vaxmap, stress_bump=transitions[:gsc_q_to_p_stress_bump])
            end
        end
    end
     

    cell.receptors .*= 0.5
    cell.receptors[vaxmap[convert_idx_to_state(cell.state[:cum_state], sim)]] = 1.0

    cell.state[:cum_state_weight] = weight
    set_cell_kind!(sim.model, cell.id, cum_state_to_idx(cell.state[:cum_state]))
end

function modular_transition_dispatcher(sim, cell, in_vaxes)
    log_cell_states(sim, cell, cell.state[:cum_state], 1.0, nothing)

    if cell.state[:resting_tracker]!=0
        cell.state[:resting_tracker] -= 1
        return
    end

    cell_type = cell.state[:type]

    transition_fn = get(TRANSITION_FN_MAP, cell_type, nothing)
    if transition_fn === nothing
        println("No transition function for cell type: $cell_type")
        return
    end

    current_state_idx = cell.state[:cum_state]

    transition_fn(sim, cell, current_state_idx, 1.0, nothing)
end

const absorb_log = Vector{Dict{Symbol, Float64}}()

function absorb_vaxes_for_cell(sim, cell, in_vaxes)
    cfg = sim.cfg
    reactions = sim.irep[:reactions]
    id = cell.id

    t = sim.model.state.t - cfg[:runtime][:burnin_steps]

    while length(absorb_log) < t
        push!(absorb_log, Dict{Symbol, Float64}())
    end

    step_log = absorb_log[t]


    cell_border_pix = get(get_stat(sim.model, :border_pixels_by_cell), id, [])
    cell_pixels = get(get_stat(sim.model, :pixels_by_cell), id, [])	

    n = cfg[:vaxes] |> length
    out_vaxes = get_outer_vaxes_avg(in_vaxes, cell_border_pix)
    #println(out_vaxes)

    for vax in 1:n
        rs = reactions_for(reactions, vax)
        for r in rs
            if r[RR_ABSORB]

                # Michaelis-Menten kinetics
                k = r[RK]
                #Km = r[RW]
                Km = 1.0
                c=0.0
                for (n, v) in r[RREACT]
                    if v == vax
                        continue
                    end
                    c = out_vaxes[v]
                end
               
                uptake = k*c / (Km + c)


                l = length(cell_pixels)
            
                for (n, v) in r[RREACT]
                    if v == vax
                        continue
                    end
                    sim.vaxes[v].x[cell_pixels] .-= n*uptake/l
                    if n*uptake != 0.0
                        vax_sym = convert_idx_to_state(v, sim)
                        step_log[vax_sym] = get(step_log, vax_sym, 0.0) + n*uptake
                    end

                end
            end
        end
    end
end


function get_outer_vaxes_avg(vaxes, cellpixels)
    map(vaxes) do vax
        vax[cellpixels] |> mean
    end
end


function plot_absorb_vaxes(sim, save_dir = nothing)
    steps = length(absorb_log)
    if steps == 0
        error("absorb_log is empty")
    end

    cell_count = length(sim.cells)
    if cell_count == 0
        error("No cells in simulation")
    end

    # Collect all vaxes that appear across the log
    vaxes = Set{Symbol}()
    for step in absorb_log
        union!(vaxes, keys(step))
    end
    vaxes = sort(collect(vaxes))

    # Build per-vax timeseries
    data = Dict{Symbol, Vector{Float64}}()
    for vax in vaxes
        data[vax] = [get(step, vax, 0.0) / cell_count for step in absorb_log]
    end

    # Plot each vax individually
    for vax in vaxes
        fig = Figure(size=(700, 400))
        ax = GLMakie.Axis(fig[1, 1], xlabel="Step", ylabel="Avg uptake per cell", title=String(vax))
        lines!(ax, 1:steps, data[vax])
        
        if save_dir !== nothing
            save(joinpath(save_dir, "$(vax)_absorption.png"), fig)
        end

        display(fig)
    end
end



function local_vax_concentration(sim, cell, vaxname::Symbol)
    vax_id = sim.irep[:vaxmap][vaxname]
    field = sim.vaxes[vax_id].x
    pixels = get(get_stat(sim.model, :border_pixels_by_cell), cell.id, [])

    if isempty(pixels)
        return 0.0  
    end

    return mean(field[p] for p in pixels)
end


function log_transition(cell_id, from_state_idx, to_state_idx, sim)
    from_state = convert_idx_to_state(from_state_idx, sim)
    to_state = convert_idx_to_state(to_state_idx, sim)
    push!(transition_log, (cell_id=cell_id, from=from_state, to=to_state, step=sim.model.state.t - sim.cfg[:runtime][:burnin_steps]))
end

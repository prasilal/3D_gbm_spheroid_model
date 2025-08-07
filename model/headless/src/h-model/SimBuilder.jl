module SimBuilder

function mk_sim!(cfg, rules = Nothing)
    index_vaxes!(cfg[:vaxes])
    vaxes = prep_vaxes_pools(cfg)

    vaxmap = vax_map(cfg[:vaxes])

    internal = Dict(
        :vaxmap => vaxmap,
        :vaxes => sort(values(cfg[:vaxes]) |> collect, by = v -> v[:id]),
        :reactions => reactions_matrix(cfg[:reactions], vaxmap),
        :zg => convert_zg_to_idxs(cfg[:rule_graph][:zg], vaxmap)
    )

    params = params_cfg(
        T = get(cfg[:rule_graph][:cpm], :T, 20),
        adhesion = AdhesionParams(adhesion_params(cfg)),
        volume = VolumeParams(volume_params(cfg)...),
        perimeter = PerimeterParams(perimeter_params(cfg)...),
        chemotaxis = map(p -> ChemotaxisParams(p...), chemotaxis_params(cfg, vaxes)),
        activity = ActivityParams(Val{:geom}(), activity_params(cfg)...),
        barrier = BarrierParams(barrier_params(cfg)),
        connectivity = ConnectivityParams(connectivity_params(cfg)),
        hard_volume_range = HardVolumeRangeParams(volume_range_params(cfg)...),
        outer_action = OuterActionParams(outer_action_params(cfg)...)
    )

    rules_cfg = mk_cfg(post_mcs_listeners = [partial(rules, cfg)])
    full_cfg = merge(StatsCfg, VolumeCfg, AdhesionCfg, PerimeterCfg, ChemotaxisCfg, ActivityCfg, BarrierCfg,
                     LocalConnectivityCfg, HardVolumeRangeCfg, OuterActionCfg,
                     params, rules_cfg)
    cpm = make_preinit_cpm(cfg[:size], cfg[:seed], cfg = full_cfg, is_torus = all(identity, cfg[:torus]))

    cfg[:sim] = Sim(Cell[], vaxes, cpm, cfg, internal)
    cfg[:sim]
end

function mk_barier!(cfg) # create barrier
    if !haskey(cfg, :barrier)
        return
    end

    vaxes = cfg[:vaxes]
    model = cfg[:sim].model
    id_barrier = CELL_TYPES_OFFSET + vaxes[cfg[:barrier][:vax]][:id]
    id = make_new_cell_id!(model, id_barrier)
    cfg[:barrier][:id] = id
    # cfg.sim.model.grid[cfg[:barrier][:pos]] .= id
    for idx in cfg[:barrier][:pos]
        setpix!(model, idx, id)
    end
end

function mk_cell!(sim, state, pos, init_receptors)
    id = seed_cell_at!(sim.model, cum_state_to_idx(state[:cum_state]), pos)
    push!(sim.cells, Cell(id, init_receptors, state))
end

function init_receptors(receptors, vaxes)
    rs = zeros(Float64, vaxes |> keys |> length)
    for (r, v) in receptors
        rs[vaxes[r][:id]] = v
    end
    rs
end

function mk_cell_from_cfg!(sim, vaxes, c)
    state = convert_state_to_idxs(c[:state], sim.irep[:vaxmap])
    mk_cell!(sim, state, CartesianIndex(c[:init_pos]), init_receptors(c[:receptors], vaxes))
end

function mk_cells!(sim)
    cfg = sim.cfg

    vaxes = cfg[:vaxes]
    for c in cfg[:cells]
        mk_cell_from_cfg!(sim, vaxes, c)
        # state = convert_state_to_idxs(c[:state], sim.irep[:vaxmap])
        # mk_cell!(sim, state, CartesianIndex(c[:init_pos]), init_receptors(c[:receptors], vaxes))
    end
end

const DEFAULT_RUNTIME_PARAMS = @d(
    :burnin_steps => 20,
    :micro_steps => 1,
    :show_fields => true,
    :show_activity => true,
    :show_sim => true,
    :cell_kind_color => default_colors,
    :field_color => default_colors_field
)

function init_sim(cfg)
    sim = mk_sim!(cfg, rules)
    mk_barier!(cfg)
    mk_cells!(sim)

    rules_enabled = get(cfg, :rules_enabled, true)

    img, t = nothing, nothing

    cfg[:runtime] = merge(DEFAULT_RUNTIME_PARAMS, get(cfg,:runtime,@d()))

    cfg[:rules_enabled] = false
    if cfg[:runtime][:show_sim]
        img,t = prep_vis_sim(sim,
                             burnin_steps = cfg[:runtime][:burnin_steps],
                             show_3d= cfg[:runtime][:show_3d],
                             num_cell_kinds = CELL_TYPES_OFFSET+length(sim.vaxes),
                             show_fields = cfg[:runtime][:show_fields],
                             show_activity = cfg[:runtime][:show_activity],
                             cell_kind_color = cfg[:runtime][:cell_kind_color],
                             field_color = cfg[:runtime][:field_color])
    else
        burnin!(sim, cfg[:runtime][:burnin_steps])
    end
    cfg[:rules_enabled] = rules_enabled

    @d(
        :sim => sim,
        :img => img,
        :t => t
    )
end

end # module

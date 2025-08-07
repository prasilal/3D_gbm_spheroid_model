module Rules

function rules(cfg, model)
    if !get(cfg, :rules_enabled, true)
        return
    end

    sim = cfg[:sim]
    in_vaxes = deepcopy(sim.vaxes)
    barrier_idxs = haskey(cfg, :barrier) ? cfg[:barrier][:pos] : nothing

    for c in sim.cells
        update_cum_state(sim, c, in_vaxes)
    end

    cells = []
    for c in sim.cells
        new_cells = update_by_zg!(sim, c, sim.vaxes)
        update_r_vaxes(sim, c)
        cells = vcat(cells, new_cells)
    end
    sim.cells = cells

    reset_stat_state!(sim.model)

    for i in 1:get(cfg, :vaxes_rule_steps, 1)
        if i > 1
            in_vaxes = deepcopy(sim.vaxes)
        end

        for r in sim.irep[:reactions]
            perform_reaction!(sim.vaxes, in_vaxes, r, cfg)
        end

        vaxes = sim.irep[:vaxes]
        for (idx,vinfo) in enumerate(vaxes)
            field = sim.vaxes[idx]
            perform_diffusion!(field, vinfo[:D], barrier_idxs)
            perform_decay!(field, vinfo[:d])
            correct_vals!(field)
        end
    end
end

end # module

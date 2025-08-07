module SimControl

using SimTypes
using Rules

always_true(_) = true

function time_macro_step!(sim)
    rules_enabled = get(sim.cfg, :rules_enabled, true)
    steps = get(sim.cfg[:runtime], :micro_steps, 1)

    sim.cfg[:rules_enabled] = false
    for _ in 1:(steps-1)
        time_step!(sim.model)
    end

    sim.cfg[:rules_enabled] = rules_enabled
    time_step!(sim.model)
end

function burnin!(sim, steps; should_we_continue = always_true)
    for i in 1:steps
        time_macro_step!(sim)
        if !should_we_continue(sim)
            break
        end
    end
end

function simulate(sim_desc; num_of_steps = 1, record = false, filename = "simulation.mp4", fps = 1.0/30., should_we_continue = always_true)
    if sim_desc[:sim].cfg[:runtime][:show_sim]
        if record
            record_sim(sim_desc[:sim], sim_desc[:t], sim_desc[:img],
                       filename = filename,
                       timestamps = 0.:fps:num_of_steps*fps,
                       should_we_continue = should_we_continue)
        else
            display(sim_desc[:img])
            anim_sim(sim_desc[:sim], sim_desc[:t],
                     num_of_steps = num_of_steps,
                     fps = fps,
                     should_we_continue = should_we_continue)
        end
    else
        burnin!(sim_desc[:sim], num_of_steps, should_we_continue = should_we_continue)
    end
    return sim_desc[:sim]
end

end # module

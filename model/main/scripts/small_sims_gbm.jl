include("hooks/custom_hooks.jl")


GRID_SIZE = 18

cfg = Dict(
    :size => (GRID_SIZE, GRID_SIZE, GRID_SIZE),
    :torus => (false, false, false),
    :seed => 1234,
    :vaxes_rule_steps => 800,

    :vaxes => @d(
        :glucose => @d(
            :D => 0.01428,       
            :d => 0.0,      
            :rd => 0.0,      
            :init => (2:17, 2:17, 2:17, 0),
            :show => true    
        ),
        :oxygen => @d(
            :D => 0.16666,       
            :d => 0.0,      
            :rd => 0.0,      
            :init => (2:17, 2:17, 2:17, 0),
            :show => true    
        ),
        :gbm_p => @d(:d => 0.0, :D => 0.0, :rd => 0.01, :diffusion_multipliers => @d(:oxygen => 0.8, :glucose => 0.7)), 
        :gbm_q => @d(:d => 0.0, :D => 0.0, :rd => 0.01, :diffusion_multipliers => @d(:oxygen => 0.5, :glucose => 0.4)),  
        :gbm_n => @d(:d => 0.0, :D => 0.0, :rd => 0.01, :diffusion_multipliers => @d(:oxygen => 0.2, :glucose => 0.1)), 
        :wall => @d(:d => 0.0, :D => 0.0, :rd => 0.01) 
    ),


    :rule_graph => @d(
        :min_weight => 0.01,
        :resting_time => 30,
        :zg => @d(
            :gbm_p => [(:gbm_p, :adhesion, 100), (:gbm_p, :volume, 1000, 40), (:gbm_p, :perimeter, 10, 500), (:gbm_p, :prod_r, 0.1)], #40
            :gbm_q => [(:gbm_q, :adhesion, 100),  (:gbm_q, :volume, 1000, 35), (:gbm_q, :perimeter, 100, 300), (:gbm_q, :prod_r, 0.1)], #35
            :gbm_n => [(:gbm_n, :adhesion, 120),  (:gbm_n, :volume, 1000, 20),  (:gbm_n, :perimeter, 100, 400), (:gbm_n, :prod_r, 0.1)], #20
            :oxygen => [(:oxygen, :adhesion, 0), (:oxygen, :volume, 0, 0), (:oxygen, :perimeter, 0, 0)],
            :glucose => [(:glucose, :adhesion, 0), (:glucose, :volume, 0, 0), (:glucose, :perimeter, 0, 0)],
            :wall => [(:wall, :adhesion, 0), (:wall, :volume, 0, 0), (:wall, :perimeter, 0, 0)]
        ),
        :cpm => @d(:T => 15, :other_adhesion => 30)
    ),

    :reactions => [
        # GBM Proliferating: high oxygen + glucose
        @d(:react => [(1, :oxygen), (1, :gbm_p)], :prod => [], :k => 5.0, :w => 1.0, :r_absorb => true),
        @d(:react => [(1, :glucose), (1, :gbm_p)], :prod => [], :k => 5.0, :w => 1.0, :r_absorb => true),

        # GBM Quiescent: reduced uptake
        @d(:react => [(1, :oxygen), (1, :gbm_q)], :prod => [], :k => 0.5, :w => 1.0, :r_absorb => true),
        @d(:react => [(1, :glucose), (1, :gbm_q)], :prod => [], :k => 0.5, :w => 1.0, :r_absorb => true),

        # GBM Necrotic: very little uptake
        @d(:react => [(1, :oxygen), (1, :gbm_n)], :prod => [], :k => 0.0, :w => 1.0, :r_absorb => true),
        @d(:react => [(1, :glucose), (1, :gbm_n)], :prod => [], :k => 0.0, :w => 1.0, :r_absorb => true)
    ],

    :cell_state_update => @d(
        :transition_fn => test_transition_dispatcher
    ),


    :cells => begin

        Random.seed!(1235)  

        num_cells = 55
        min_dist = 3  
        cells = []

        cell_centers = []

        function is_far_enough(pos, others, min_dist)
            all(p -> norm(ntuple(i -> p[i] - pos[i], 3)) ≥ min_dist, others)
        end

        while length(cell_centers) < num_cells
            candidate = (rand(3:16), rand(3:16), rand(3:16))  
            if is_far_enough(candidate, cell_centers, min_dist)
                push!(cell_centers, candidate)
            end
        end

        for center in cell_centers
            push!(cells, @d(
                :state => @d(:type => :gbm, :cum_state => :gbm_p, :cum_state_weight => 1.0),
                :init_pos => center,
                :receptors => @d(:gbm_p => 1.0)
            ))
        end

        cells
    end,

    :runtime => @d(
        :show_sim => true,
        :show_3d => true,
        :show_fields => true,
        :show_activity => false,
        :burnin_steps => 50,
        :micro_steps => 1
    ),

    :outer_action => @d(
        :regions => [((17, 2, 2), (17, 17, 17), :inject), ((2, 2, 2), (2, 17, 17), :absorb) ],
        :injection_amount => @d(:glucose => 4.33, :oxygen => 5.94),
    ),

    :barrier => @d(
        :vax => :wall,  
        :pos => vcat(
            vec([CartesianIndex(x, y, 1)     for x in 1:18, y in 1:18]),
            vec([CartesianIndex(x, y, 18)    for x in 1:18, y in 1:18]),
            vec([CartesianIndex(x, 1, z)     for x in 1:18, z in 1:18]),
            vec([CartesianIndex(x, 18, z)    for x in 1:18, z in 1:18]),
            vec([CartesianIndex(1, y, z)     for y in 1:18, z in 1:18]),
            vec([CartesianIndex(18, y, z)    for y in 1:18, z in 1:18])
        )
    )
    
    
)
cfg[:outer_action][:action] = partial(inject_absorb_track_function, cfg)

sim_desc = init_sim(cfg)


simulate(sim_desc, record= false, num_of_steps = 5)

glucose_field = cfg[:sim].vaxes[cfg[:vaxes][:glucose][:id]].x
oxygen_field = cfg[:sim].vaxes[cfg[:vaxes][:oxygen][:id]].x

#plot_absorb_vaxes(cfg[:sim], "outputs")


# Plot the outflow history
plot_outflow_history(:glucose, "outputs/glucose_outflow.png")
plot_outflow_history(:oxygen, "outputs/oxygen_outflow.png")

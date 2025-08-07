include("../hooks/custom_hooks.jl")


cfg = Dict(
    :size => (42, 42, 42),
    :torus => (false, false, false),
    :seed => 1234,
    :vaxes_rule_steps => 500,

    :vaxes => @d(
        :glucose => @d(
            :D => 0.1667,       #  Effective diffusion
            :d => 0.0,      # Decay
            :rd => 0.0,      
            :init => (2:41, 2:41, 2:41, 0.0),
            :show => true    
        ),
        :oxygen => @d(
            :D => 0.8,       # Diffusion
            :d => 0.0,      # Decay
            :rd => 0.0,      
            :init => (2:41, 2:41, 2:41, 0.0),
            :show => true    
        ),
        :gbm_p => @d(:d => 0.0, :D => 0.0, :rd => 0.01, :diffusion_multipliers => @d(:oxygen => 0.6, :glucose => 0.5)), #0.1
        :gbm_q => @d(:d => 0.0, :D => 0.0, :rd => 0.01, :diffusion_multipliers => @d(:oxygen => 0.4, :glucose => 0.3)), #0.2 
        :gbm_n => @d(:d => 0.0, :D => 0.0, :rd => 0.01, :diffusion_multipliers => @d(:oxygen => 0.2, :glucose => 0.1)), #0.05
        :wall => @d(:d => 0.0, :D => 0.0, :rd => 0.01) # Wall
    ),


    :rule_graph => @d(
        :min_weight => 0.01,
        :resting_time => 30,
        :zg => @d(
            :gbm_p => [(:gbm_p, :adhesion, 20), (:gbm_p, :volume, 100, 144), (:gbm_p, :perimeter, 100, 600), (:gbm_p, :prod_r, 0.1), (:gbm_p, :divide, 5)],
            :gbm_q => [(:gbm_q, :adhesion, 40),  (:gbm_q, :volume, 100, 80), (:gbm_q, :perimeter, 100, 470), (:gbm_q, :prod_r, 0.1)],
            :gbm_n => [(:gbm_n, :adhesion, 50),  (:gbm_n, :volume, 100, 65),  (:gbm_n, :perimeter, 100, 420), (:gbm_n, :prod_r, 0.1)],
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
        @d(:react => [(1, :oxygen), (1, :gbm_q)], :prod => [], :k => 0.0003, :w => 0.5, :r_absorb => true),
        @d(:react => [(1, :glucose), (1, :gbm_q)], :prod => [], :k => 0.0006, :w => 0.5, :r_absorb => true),

        # GBM Necrotic: very little uptake
        @d(:react => [(1, :oxygen), (1, :gbm_n)], :prod => [], :k => 1.0, :w => 1.0, :r_absorb => true),
        @d(:react => [(1, :glucose), (1, :gbm_n)], :prod => [], :k => 1.0, :w => 1.0, :r_absorb => true)
    ],

    :cell_state_update => @d(
        :transition_fn => test_transition_dispatcher
    ),


    :cells => begin

        Random.seed!(1235)  

        num_cells = 280
        min_dist = 4  
        cells = []

        cell_centers = []

        function is_far_enough(pos, others, min_dist)
            all(p -> norm(ntuple(i -> p[i] - pos[i], 3)) ≥ min_dist, others)
        end

        while length(cell_centers) < num_cells
            candidate = (rand(3:40), rand(3:40), rand(3:40))  
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
        :regions => [((41, 2, 2), (41, 41, 41), :inject), ((2, 2, 2), (2, 41, 41), :absorb) ],
        :injection_amount => @d(:glucose => 0.5, :oxygen => 0.5),
    ),

    :barrier => @d(
        :vax => :wall,  
        :pos => vcat(
            vec([CartesianIndex(x, y, 1)     for x in 1:42, y in 1:42]),
            vec([CartesianIndex(x, y, 42)    for x in 1:42, y in 1:42]),
            vec([CartesianIndex(x, 1, z)     for x in 1:42, z in 1:42]),
            vec([CartesianIndex(x, 42, z)    for x in 1:42, z in 1:42]),
            vec([CartesianIndex(1, y, z)     for y in 1:42, z in 1:42]),
            vec([CartesianIndex(42, y, z)    for y in 1:42, z in 1:42])
        )
    )
    
    
)
cfg[:outer_action][:action] = partial(inject_absorb_track_function, cfg)

sim_desc = init_sim(cfg)


simulate(sim_desc, record= false, num_of_steps = 30)

glucose_field = cfg[:sim].vaxes[cfg[:vaxes][:glucose][:id]].x
oxygen_field = cfg[:sim].vaxes[cfg[:vaxes][:oxygen][:id]].x
#plot_line_profile(glucose_field, 21, 21)
#plot_averaged_x_profile(glucose_field)
plot_field_slice_xy(glucose_field, 21)
plot_field_slice_xy(oxygen_field, 21)


"""
# Plot the outflow history
plot_outflow_history(:glucose)
plot_outflow_history(:oxygen)"""
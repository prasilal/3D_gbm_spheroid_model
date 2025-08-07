include("../hooks/custom_hooks.jl")


cfg = Dict(
    :size => (30, 30, 30),
    :torus => (false, false, false),
    :seed => 1234,

    :vaxes => @d(
        :glucose => @d(
            :D => 0.5,       # Diffusion
            :d => 0.01,      # Decay
            :rd => 0.0,      # Receptor decay (not used here)
            :init => (:, :, 1:15, 1.0),
            :show => true,    # Show in visualization
            :thresholds => [0.01,0.065]
        ),
        :oxygen => @d(
            :D => 0.5,       # Diffusion
            :d => 0.01,      # Decay
            :rd => 0.0,      # Receptor decay (not used here)
            :init => (:, :, 1:15, 1.0),
            :show => true,    # Show in visualization
            :thresholds => [0.01,0.065]
        ),
        :gbm_p => @d(:d => 0.0, :D => 0.0, :rd => 0.01, :diffusion_multipliers => @d(:oxygen => 0.05, :glucose => 0.0), :interphase_duration => 20), #0.1
        :gbm_q => @d(:d => 0.0, :D => 0.0, :rd => 0.01, :diffusion_multipliers => @d(:oxygen => 0.1, :glucose => 0.0), :interphase_duration => 50), #0.2 
        :gbm_n => @d(:d => 0.0, :D => 0.0, :rd => 0.01, :diffusion_multipliers => @d(:oxygen => 0.01, :glucose => 0.0), :interphase_duration => 1000) #0.05
    ),


    :rule_graph => @d(
        :min_weight => 0.01,
        :resting_time => 5,
        :zg => @d(
            :gbm_p => [(:gbm_p, :adhesion, 1000), (:gbm_p, :volume, 100, 15), (:gbm_p, :perimeter, 100, 250), (:gbm_p, :prod_r, 0.01), (:gbm_p, :divide, 0.9)],
            :gbm_q => [(:gbm_q, :adhesion, 800),  (:gbm_q, :volume, 100, 13), (:gbm_q, :perimeter, 100, 220), (:gbm_q, :prod_r, 0.01), (:gbm_q, :divide, 0.5)],
            :gbm_n => [(:gbm_n, :adhesion, 2000),  (:gbm_n, :volume, 50, 10),  (:gbm_n, :perimeter, 80, 150)],
            :oxygen => [(:oxygen, :adhesion, 0), (:oxygen, :volume, 0, 0), (:oxygen, :perimeter, 0, 0)],
            :glucose => [(:glucose, :adhesion, 0), (:glucose, :volume, 0, 0), (:glucose, :perimeter, 0, 0)]
        ),
        :cpm => @d(:T => 15, :other_adhesion => 1000),
        :division_hook => test_division_hook
    ),

    :reactions => [
        # GBM Proliferating: high oxygen + glucose
        @d(:react => [(1, :oxygen), (1, :gbm_p)], :prod => [(1, :gbm_p)], :k => 0.002, :w => 1.0, :r_absorb => true),
        @d(:react => [(1, :glucose), (1, :gbm_p)], :prod => [(1, :gbm_p)], :k => 0.003, :w => 1.0, :r_absorb => true),

        # GBM Quiescent: reduced uptake
        @d(:react => [(1, :oxygen), (1, :gbm_q)], :prod => [], :k => 0.0003, :w => 0.5, :r_absorb => true),
        @d(:react => [(1, :glucose), (1, :gbm_q)], :prod => [], :k => 0.0006, :w => 0.5, :r_absorb => true),

        # GBM Necrotic: very little uptake
        @d(:react => [(1, :oxygen), (1, :gbm_n)], :prod => [], :k => 0.00001, :w => 0.1, :r_absorb => true),
        @d(:react => [(1, :glucose), (1, :gbm_n)], :prod => [], :k => 0.00002, :w => 0.1, :r_absorb => true)
    ],

    :cell_state_update => @d(
        :transition_fn => test_transition_dispatcher
    ),


    :cells => begin
        # Create a 3D spheroid-like cluster at the center
        spheroid_center = (15, 15, 15)
        spheroid_radius = 4
        cells = []
        
        for x in (spheroid_center[1]-spheroid_radius):(spheroid_center[1]+spheroid_radius), 
            y in (spheroid_center[2]-spheroid_radius):(spheroid_center[2]+spheroid_radius),
            z in (spheroid_center[3]-spheroid_radius):(spheroid_center[3]+spheroid_radius)
        
            if sqrt((x - spheroid_center[1])^2 + (y - spheroid_center[2])^2 + (z - spheroid_center[3])^2) <= spheroid_radius
                push!(cells, @d(
                    :state => @d(:type => :gbm, :cum_state => :gbm_p, :cum_state_weight => 1.0, :phase => :interphase, :cell_cycle_timer => 1),
                    :init_pos => (x, y, z),
                    :receptors => @d(:gbm_p => 1.0)
                ))
            end
        end
        cells
    end,

    :runtime => @d(
        :show_sim => true,
        :show_3d => true,
        :show_fields => true,
        :show_activity => false,
        :burnin_steps => 20
    )
)

sim_desc = init_sim(cfg)
simulate(sim_desc, record= false, num_of_steps = 10)

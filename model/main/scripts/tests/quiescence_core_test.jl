include("../hooks/custom_hooks.jl")



cfg = Dict(
    :size => (30, 30, 30),
    :torus => (false, false, false),
    :seed => 1234,

    :vaxes => @d(
        :glucose => @d(
            :D => 0.1,       # Diffusion
            :d => 0.01,      # Decay
            :rd => 0.0,      # Receptor decay (not used here)
            :init => (1:15,:,:, 0.7),
            :show => true    # Show in visualization
        ),
        :oxygen => @d(
            :D => 0.1,       # Diffusion
            :d => 0.01,      # Decay
            :rd => 0.0,      # Receptor decay (not used here)
            :init => (1:15,:,:, 0.7),
            :show => true    # Show in visualization
        ),
        :gbm_p => @d(:d => 0.0, :D => 0.0, :rd => 0.0),
        :gbm_q => @d(:d => 0.0, :D => 0.0, :rd => 0.0),
        :gbm_n => @d(:d => 0.0, :D => 0.0, :rd => 0.0)
    ),


    :rule_graph => @d(
        :min_weight => 0.01,
        :resting_time => 5,
        :zg => @d(
            :gbm_p => [(:gbm_p, :adhesion, 150), (:gbm_p, :volume, 100, 150), (:gbm_p, :perimeter, 100, 250), (:gbm_p, :prod_r, 0.01), (:gbm_p, :divide, 5, (:gbm_p, 80))],
            :gbm_q => [(:gbm_q, :adhesion, 150),  (:gbm_q, :volume, 100, 130), (:gbm_q, :perimeter, 100, 220), (:gbm_q, :prod_r, 0.005)],
            :gbm_n => [(:gbm_n, :adhesion, 150),  (:gbm_n, :volume, 50, 100),  (:gbm_n, :perimeter, 80, 150)],
            :oxygen => [(:oxygen, :adhesion, 0), (:oxygen, :volume, 0, 0), (:oxygen, :perimeter, 0, 0)],
            :glucose => [(:glucose, :adhesion, 0), (:glucose, :volume, 0, 0), (:glucose, :perimeter, 0, 0)]
        ),
        :cpm => @d(:T => 15, :other_adhesion => 30) 
    ),

    :reactions => [
        @d(:react => [(1, :oxygen), (1, :glucose)], :prod => [(1, :gbm_p)], :k => 0.8, :w => 1.0, :r_absorb => true)
    ],

    :cell_state_update => @d(
        :transition_fn => test_transition_dispatcher
    ),


    :cells => begin
        spheroid_center = (15, 15, 15)
        spheroid_radius = 4
        cells = []
        
        for x in (spheroid_center[1]-spheroid_radius):(spheroid_center[1]+spheroid_radius), 
            y in (spheroid_center[2]-spheroid_radius):(spheroid_center[2]+spheroid_radius),
            z in (spheroid_center[3]-spheroid_radius):(spheroid_center[3]+spheroid_radius)
        
            if sqrt((x - spheroid_center[1])^2 + (y - spheroid_center[2])^2 + (z - spheroid_center[3])^2) <= spheroid_radius
                push!(cells, @d(
                    :state => @d(:type => :gbm, :cum_state => :gbm_p, :cum_state_weight => 1.0),
                    :init_pos => (x, y, z),
                    :receptors => @d(:glucose => 0.0)
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
        :burnin_steps => 5
    )
)

sim_desc = init_sim(cfg)
simulate(sim_desc, num_of_steps = 10, fps = 0.5)
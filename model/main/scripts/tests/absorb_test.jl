include("../hooks/custom_hooks.jl")

grid_size_test = 10

cfg = Dict(
    :size => (grid_size_test, grid_size_test, grid_size_test),
    :torus => (false, false, false),
    :seed => 1234,
    :vaxes_rule_steps => 50,

    :vaxes => @d(
        :glucose => @d(
            :D => 0.1667,       #  Effective diffusion
            :d => 0.0,      # Decay
            :rd => 0.0,      
            :init => (2:grid_size_test-1, 2:grid_size_test-1, 2:grid_size_test-1, 0.8),
            :show => true    
        ),
        :oxygen => @d(
            :D => 0.1667,       # Diffusion
            :d => 0.0,      # Decay
            :rd => 0.0,      
            :init => (2:grid_size_test-1, 2:grid_size_test-1, 2:grid_size_test-1, 0.8),
            :show => true    
        ),
        :gbm_p => @d(:d => 0.0, :D => 0.0, :rd => 0.01),
        :wall => @d(:d => 0.0, :D => 0.0, :rd => 0.01) # Wall
    ),


    :rule_graph => @d(
        :min_weight => 0.01,
        :resting_time => 30,
        :zg => @d(
            :gbm_p => [(:gbm_p, :adhesion, 40), (:gbm_p, :volume, 100, 20), (:gbm_p, :perimeter, 100, 100), (:gbm_p, :prod_r, 0.1), (:gbm_p, :divide, 5)],
            :oxygen => [(:oxygen, :adhesion, 0), (:oxygen, :volume, 0, 0), (:oxygen, :perimeter, 0, 0)],
            :glucose => [(:glucose, :adhesion, 0), (:glucose, :volume, 0, 0), (:glucose, :perimeter, 0, 0)],
            :wall => [(:wall, :adhesion, 0), (:wall, :volume, 0, 0), (:wall, :perimeter, 0, 0)]
        ),
        :cpm => @d(:T => 15, :other_adhesion => 30)
    ),

    :reactions => [
        # GBM Proliferating: high oxygen + glucose
        @d(:react => [(1, :oxygen), (1, :gbm_p)], :prod => [], :k => 1.0, :w => 1.0, :r_absorb => true),
        @d(:react => [(1, :glucose), (1, :gbm_p)], :prod => [], :k => 0.9, :w => 1.0, :r_absorb => true)
    ],

    :cell_state_update => @d(
        :transition_fn => test_transition_dispatcher
    ),


    :cells => begin

        Random.seed!(1235)  

        num_cells = 3
        min_dist = 4  
        cells = []

        cell_centers = []

        function is_far_enough(pos, others, min_dist)
            all(p -> norm(ntuple(i -> p[i] - pos[i], 3)) ≥ min_dist, others)
        end

        while length(cell_centers) < num_cells
            candidate = (rand(2:grid_size_test-1), rand(2:grid_size_test-1), rand(2:grid_size_test-1))  
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

    :barrier => @d(
        :vax => :wall,  
        :pos => vcat(
            vec([CartesianIndex(x, y, 1)     for x in 1:grid_size_test, y in 1:grid_size_test]),
            vec([CartesianIndex(x, y, grid_size_test)    for x in 1:grid_size_test, y in 1:grid_size_test]),
            vec([CartesianIndex(x, 1, z)     for x in 1:grid_size_test, z in 1:grid_size_test]),
            vec([CartesianIndex(x, grid_size_test, z)    for x in 1:grid_size_test, z in 1:grid_size_test]),
            vec([CartesianIndex(1, y, z)     for y in 1:grid_size_test, z in 1:grid_size_test]),
            vec([CartesianIndex(grid_size_test, y, z)    for y in 1:grid_size_test, z in 1:grid_size_test])
        )
    )

    
    
)

sim_desc = init_sim(cfg)


simulate(sim_desc, record= false, num_of_steps = 50)

plot_absorb_vaxes(cfg[:sim])

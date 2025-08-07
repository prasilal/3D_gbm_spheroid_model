include("hooks/custom_hooks.jl")

function mk_cfg(args)
    GRID_SIZE = args["gridsize"]
    NUM_CELLS = args["nofcells"]
    SEED = args["seed"]
    GBM_RATIO = args["gbm"]
    GSC_RATIO = args["gsc"]
    MSC_RATIO = args["msc"]

    cfg = Dict(
        :size => (GRID_SIZE, GRID_SIZE, GRID_SIZE),
        :torus => (false, false, false),
        :seed => SEED,
        :vaxes_rule_steps => 2093,

        :vaxes => @d(
            :glucose => @d(
                :D => 0.01428,       
                :d => 0.0,      
                :rd => 0.0,      
                :init => (2:GRID_SIZE-1, 2:GRID_SIZE-1, 2:GRID_SIZE-1, 1.6037),
                :show => false    
            ),
            :oxygen => @d(
                :D => 0.16666,       
                :d => 0.0,      
                :rd => 0.0,      
                :init => (2:GRID_SIZE-1, 2:GRID_SIZE-1, 2:GRID_SIZE-1, 2.2),
                :show => false    
            ),
            :lactate => @d(
                :D => 0.005,       
                :d => 0.00001,      
                :rd => 0.0,      
                :init => (2:GRID_SIZE-1, 2:GRID_SIZE-1, 2:GRID_SIZE-1, 0),
                :show => true    
            ),
            :ecm_gsc => @d(
                :D => 0.0001,       
                :d => 0.00001,      
                :rd => 0.0,      
                :init => (2:GRID_SIZE-1, 2:GRID_SIZE-1, 2:GRID_SIZE-1, 0),
                :show => false    
            ),
            :ecm_barrier => @d(
                :D => 0.0001,       
                :d => 0.00001,      
                :rd => 0.0,      
                :init => (2:GRID_SIZE-1, 2:GRID_SIZE-1, 2:GRID_SIZE-1, 0),
                :show => false    
            ),
            :growth_signal => @d(
                :D => 0.001,       
                :d => 0.0001,      
                :rd => 0.0,      
                :init => (2:GRID_SIZE-1, 2:GRID_SIZE-1, 2:GRID_SIZE-1, 0),
                :show => false    
            ),
            :gbm_p => @d(:d => 0.0, :D => 0.0, :rd => 0.0, :diffusion_multipliers => @d(:oxygen => 0.8, :glucose => 0.7), :interphase_duration => 10, :resting_period => 20), 
            :gbm_q => @d(:d => 0.0, :D => 0.0, :rd => 0.0, :diffusion_multipliers => @d(:oxygen => 0.5, :glucose => 0.4), :interphase_duration => 15, :resting_period => 20),  
            :gbm_n => @d(:d => 0.0, :D => 0.0, :rd => 0.0, :diffusion_multipliers => @d(:oxygen => 0.2, :glucose => 0.1)),

            :gsc_p => @d(:d => 0.0, :D => 0.0, :rd => 0.0, :diffusion_multipliers => @d(:oxygen => 0.8, :glucose => 0.7), :interphase_duration => 10, :resting_period => 30),
            :gsc_q => @d(:d => 0.0, :D => 0.0, :rd => 0.0, :diffusion_multipliers => @d(:oxygen => 0.7, :glucose => 0.5), :interphase_duration => 10, :resting_period => 30),

            :msc_p => @d(:d => 0.0, :D => 0.0, :rd => 0.0, :diffusion_multipliers => @d(:oxygen => 0.7, :glucose => 0.5), :interphase_duration => 15, :resting_period => 100), 
            :msc_anti => @d(:d => 0.0, :D => 0.0, :rd => 0.0, :diffusion_multipliers => @d(:oxygen => 0.7, :glucose => 0.5), :interphase_duration => 15, :resting_period => 15),
            :msc_before => @d(:d => 0.0, :D => 0.0, :rd => 0.0, :diffusion_multipliers => @d(:oxygen => 0.5, :glucose => 0.4), :interphase_duration => 20, :resting_period => 15),
            :wall => @d(:d => 0.0, :D => 0.0, :rd => 0.01) 
        ),


        :rule_graph => @d(
            :min_weight => 0.01,
            :resting_time => 30,
            :zg => @d(
                :gbm_p => [(:gbm_p, :adhesion, 100), (:gbm_p, :volume, 2000, 11), (:gbm_p, :perimeter, 10, 140), (:gbm_p, :prod_r, 0.1), (:lactate, :prod_v, 1.0), (:gbm_p, :divide, 0.7)], #40
                :gbm_q => [(:gbm_q, :adhesion, 100),  (:gbm_q, :volume, 1000, 10), (:gbm_q, :perimeter, 10, 130), (:gbm_q, :prod_r, 0.1), (:gbm_q, :divide, 0.2)], #35
                :gbm_n => [(:gbm_n, :adhesion, 120),  (:gbm_n, :volume, 1000, 10),  (:gbm_n, :perimeter, 10, 130), (:gbm_n, :prod_r, 0.1)], #20

                :gsc_p => [(:gsc_p, :adhesion, 80), (:gsc_p, :volume, 2000, 10), (:gsc_p, :perimeter, 10, 125), (:gsc_p, :prod_r, 0.1), (:lactate, :prod_v, 1.0), (:gsc_p, :divide, 0.75)],
                :gsc_q => [(:gsc_q, :adhesion, 60), (:gsc_q, :volume, 1000, 10), (:gsc_q, :perimeter, 10, 125), (:gsc_q, :prod_r, 0.1)],

                :msc_p => [(:msc_p, :adhesion, 1000), (:msc_p, :volume, 2000, 10), (:msc_p, :perimeter, 10, 120), (:msc_p, :prod_r, 0.1), (:ecm_gsc, :prod_v, 1.5), (:growth_signal, :prod_v, 1.5)],
                :msc_anti => [(:msc_anti, :adhesion, 1000), (:msc_anti, :volume, 2000, 10), (:msc_anti, :perimeter, 10, 120), (:msc_anti, :prod_r, 0.1), (:ecm_gsc, :prod_v, 1.5), (:growth_signal, :prod_v, 1.5)],
                :msc_before => [(:msc_before, :adhesion, 600), (:msc_before, :volume, 1000, 15), (:msc_before, :perimeter, 10, 230), (:msc_before, :prod_r, 0.1)],

                :lactate => [(:lactate, :adhesion, 0), (:lactate, :volume, 0, 0), (:lactate, :perimeter, 0, 0)],
                :ecm_gsc => [(:ecm_gsc, :adhesion, 0), (:ecm_gsc, :volume, 0, 0), (:ecm_gsc, :perimeter, 0, 0)],
                :ecm_barrier => [(:ecm_barrier, :adhesion, 0), (:ecm_barrier, :volume, 0, 0), (:ecm_barrier, :perimeter, 0, 0)],
                :growth_signal => [(:growth_signal, :adhesion, 0), (:growth_signal, :volume, 0, 0), (:growth_signal, :perimeter, 0, 0)],
                :oxygen => [(:oxygen, :adhesion, 0), (:oxygen, :volume, 0, 0), (:oxygen, :perimeter, 0, 0)],
                :glucose => [(:glucose, :adhesion, 0), (:glucose, :volume, 0, 0), (:glucose, :perimeter, 0, 0)],
                :wall => [(:wall, :adhesion, 0), (:wall, :volume, 0, 0), (:wall, :perimeter, 0, 0)]
            ),
            :cpm => @d(:T => 15, :other_adhesion => 30),
            :division_hook => division_hook
        ),

        :reactions => [
            # GBM Proliferating: high oxygen + glucose
            @d(:react => [(1, :oxygen), (1, :gbm_p)], :prod => [], :k => 5.0, :w => 1.0, :r_absorb => true),
            @d(:react => [(1, :glucose), (1, :gbm_p)], :prod => [], :k => 5.0, :w => 1.0, :r_absorb => true),
            # GBM Quiescent: reduced uptake
            @d(:react => [(1, :oxygen), (1, :gbm_q)], :prod => [], :k => 0.5, :w => 1.0, :r_absorb => true),
            @d(:react => [(1, :glucose), (1, :gbm_q)], :prod => [], :k => 0.5, :w => 1.0, :r_absorb => true),
            # GBM Necrotic: no uptake
            @d(:react => [(1, :oxygen), (1, :gbm_n)], :prod => [], :k => 0.0, :w => 1.0, :r_absorb => true),
            @d(:react => [(1, :glucose), (1, :gbm_n)], :prod => [], :k => 0.0, :w => 1.0, :r_absorb => true),

            # GSC Proliferating: high oxygen + glucose
            @d(:react => [(1, :oxygen), (1, :gsc_p)], :prod => [], :k => 5.0, :w => 1.0, :r_absorb => true),
            @d(:react => [(1, :glucose), (1, :gsc_p)], :prod => [], :k => 5.0, :w => 1.0, :r_absorb => true),
            # GSC Quiescent: reduced uptake 
            @d(:react => [(1, :oxygen), (1, :gsc_q)], :prod => [], :k => 0.5, :w => 1.0, :r_absorb => true),
            @d(:react => [(1, :glucose), (1, :gsc_q)], :prod => [], :k => 0.5, :w => 1.0, :r_absorb => true),

            # MSC Pro-tumor: 
            @d(:react => [(1, :oxygen), (1, :msc_p)], :prod => [], :k => 1.0, :w => 1.0, :r_absorb => true),
            @d(:react => [(1, :glucose), (1, :msc_p)], :prod => [], :k => 1.0, :w => 1.0, :r_absorb => true),
            # MSC Anti-tumor: 
            @d(:react => [(1, :oxygen), (1, :msc_anti)], :prod => [], :k => 1.0, :w => 1.0, :r_absorb => true),
            @d(:react => [(1, :glucose), (1, :msc_anti)], :prod => [], :k => 1.0, :w => 1.0, :r_absorb => true),
            # MSC Before:
            @d(:react => [(1, :oxygen), (1, :msc_before)], :prod => [], :k => 0.7, :w => 1.0, :r_absorb => true),
            @d(:react => [(1, :glucose), (1, :msc_before)], :prod => [], :k => 0.7, :w => 1.0, :r_absorb => true)
        ],

        :cell_state_update => @d(
            :transition_fn => modular_transition_dispatcher,
            :absorb_fn => absorb_vaxes_for_cell,
            :thresholds => @d(
                :hypoxia => 0.35,
                :severe_hypoxia => 0.05,
                :starvation => 0.2962,
                :severe_starvation => 0.074,
                :growth_signal => 3.0,
                :ecm_gsc => 4.0,
                :env_score_threshold => @d(:low => 0.2, :medium => 0.3, :high => 0.4),
                :acidity => 4.0
            ),
            :transitions => @d(

                # GBM thresholds
                :gbm_stress_q => 5,
                :gbm_stress_n => 10,
                :gbm_p_to_n_prob => 0.95,
                :gbm_p_to_q_prob => @d(:verylow => 0.65, :low => 0.75, :medium => 0.85, :high => 0.95),
                :gbm_to_gsc_prob => 0.80,
                :gbm_q_to_n_prob => @d(:low => 0.70, :medium => 0.80, :high => 0.95),
                :gbm_q_to_p_prob => @d(:low => 0.30, :medium => 0.40, :high => 0.90),
                :gbm_q_to_p_stress_bump => 3,
            
                # GSC thresholds
                :gsc_stress_q => 10,
                :gsc_stress_n => 15,
                :gsc_p_to_gbm_prob => 0.95,
                :gsc_p_to_q_prob => @d(:verylow => 0.65, :low => 0.70, :medium => 0.75, :high => 0.80, :higher => 0.90, :veryhigh => 0.95),
                :gsc_q_to_gbm_prob => @d(:low => 0.60, :medium => 0.70, :high => 0.85),
                :gsc_q_to_p_prob => @d(:low => 0.30, :medium => 0.40, :high => 0.90),
                :gsc_q_to_p_stress_bump => 5,
            
                # MSC thresholds
                :msc_stress_q => 5,
                :msc_stress_n => 10,
                :msc_p_to_anti_prob => @d(:verylow => 0.4, :low => 0.5, :medium => 0.74, :high => 0.90),
                :msc_anti_to_p_prob => @d(:low => 0.5, :medium => 0.75, :high => 0.95),
                :gbm_contact_prob => @d(:verylow => 0.1, :low => 0.3, :medium => 0.5, :high => 0.7),
                :msc_p_to_anti_stress_bump => 2,
                :acidity_thresholds => @d(:low => 0.3, :medium => 0.5, :high => 0.7)
            )
        ),


        :cells => begin
            Random.seed!(SEED)

            num_cells = NUM_CELLS
            ratios = Dict(:gbm => GBM_RATIO, :gsc => GSC_RATIO, :msc => MSC_RATIO)
            min_dist = 1
            local i = 1

            # Estimate center and spheroid radius
            local center = (GRID_SIZE÷2, GRID_SIZE÷2, GRID_SIZE÷2)
            if MSC_RATIO == 1.0
                r_spheroid = 18
            else
                r_spheroid = 20
            end

            # Assign how many cells per type
            counts = Dict(k => round(Int, ratios[k] * num_cells) for k in keys(ratios))
            total = sum(values(counts))

            # Prepare cell centers
            cell_centers = []
            function is_far_enough(pos, others, min_dist)
                all(p -> norm(ntuple(i -> p[i] - pos[i], 3)) ≥ min_dist, others)
            end
            while length(cell_centers) < total
                candidate = ntuple(i -> rand(center[i] - r_spheroid:center[i] + r_spheroid), 3)
                if norm(SVector(candidate...) .- SVector(center...)) ≤ r_spheroid &&
                is_far_enough(candidate, cell_centers, min_dist)
                    push!(cell_centers, candidate)
                end
            end

            # Assign types to cells in order
            cells = []
            for (ctype, count) in counts
                for _ in 1:count
                    center = cell_centers[i]
                    cum_state = ctype == :msc ? :msc_before : Symbol(string(ctype), "_p")
                    push!(cells, @d(
                        :state => @d(:type => ctype, :cum_state => cum_state, :cum_state_weight => 1.0, :phase => :interphase, :cell_cycle_timer => rand(1:10), :stress => 0, :resting_tracker => rand(5:30)),
                        :init_pos => center,
                        :receptors => @d(cum_state => 1.0)
                    ))
                
                    i+=1
                end
            end
            

            cells
        end,


        :runtime => @d(
            :show_sim => false,
            :show_3d => false,
            :show_fields => false,
            :show_activity => false,
            :burnin_steps => 20,
            :micro_steps => 1
        ),

        :outer_action => @d(
            :regions => [((GRID_SIZE-1, 2, 2), (GRID_SIZE-1, GRID_SIZE-1, GRID_SIZE-1), :inject), 
            ((2, 2, 2), (GRID_SIZE-1, 2, GRID_SIZE-1), :inject),
            ((2, 2, 2), (2, GRID_SIZE-1, GRID_SIZE-1), :inject),
            ((2, GRID_SIZE-1, 2), (GRID_SIZE-1, GRID_SIZE-1, GRID_SIZE-1), :inject), 
            ((2, 2, 2), (GRID_SIZE-1, GRID_SIZE-1, 2), :inject), 
            ((2, 2, GRID_SIZE-1), (GRID_SIZE-1, GRID_SIZE-1, GRID_SIZE-1), :inject)],
            :injection_amount => @d(:glucose => 1.6037, :oxygen => 2.2),
        ),

        :barrier => @d(
            :vax => :wall,  
            :pos => vcat(
                vec([CartesianIndex(x, y, 1)     for x in 1:GRID_SIZE, y in 1:GRID_SIZE]),
                vec([CartesianIndex(x, y, GRID_SIZE)    for x in 1:GRID_SIZE, y in 1:GRID_SIZE]),
                vec([CartesianIndex(x, 1, z)     for x in 1:GRID_SIZE, z in 1:GRID_SIZE]),
                vec([CartesianIndex(x, GRID_SIZE, z)    for x in 1:GRID_SIZE, z in 1:GRID_SIZE]),
                vec([CartesianIndex(1, y, z)     for y in 1:GRID_SIZE, z in 1:GRID_SIZE]),
                vec([CartesianIndex(GRID_SIZE, y, z)    for y in 1:GRID_SIZE, z in 1:GRID_SIZE])
            )
        )
        
        
    )
    cfg[:outer_action][:action] = partial(inject_absorb_track_function, cfg)

    return cfg
end

function run_sim(args)
    cfg = mk_cfg(args)
    sim_desc = init_sim(cfg)

    simulate(sim_desc, record= false, num_of_steps = args["nofsteps"])
end

module CPMParams

function adhesion_params(cfg)
    cpm = cfg[:rule_graph][:cpm]
    vaxes = cfg[:vaxes]
    if haskey(cpm, :adhesion)
        return cpm[:adhesion]
    end

    n = types_cnt(cfg)
    adh_mtx = fill(convert(Float64, get(cpm, :other_adhesion, 0.)), (n,n))
    adh_mtx[1,1] = 0.0

    zg = cfg[:rule_graph][:zg]
    for v in zg |> keys
        idx = findfirst(x -> x[2] == :adhesion, zg[v])
        if idx != nothing
            i = CELL_TYPES_OFFSET + vaxes[v][:id]
            adh_mtx[i, i] = zg[v][idx][3]
        end
    end

    for i in 2:n
        for j in 2:n
            if i != j
                adh_mtx[i,j] = min(adh_mtx[i,i], adh_mtx[j,j])
            end
        end
    end

    adh_mtx
end

function cpm_params(cfg, key)
    cpm = cfg[:rule_graph][:cpm]
    vaxes = cfg[:vaxes]
    if haskey(cpm, key)
        return cpm[key]
    end

    n = types_cnt(cfg)
    w = zeros(Float64, n)
    val = zeros(Float64, n)

    zg = cfg[:rule_graph][:zg]
    for v in zg |> keys
        idx = findfirst(x -> x[2] == key, zg[v])
        if idx != nothing
            i = CELL_TYPES_OFFSET + vaxes[v][:id]
            w[i] = zg[v][idx][3]
            val[i] = zg[v][idx][4]
        end
    end

    (w, val)
end

volume_params(cfg) = cpm_params(cfg, :volume)
perimeter_params(cfg) = cpm_params(cfg, :perimeter)
activity_params(cfg) = cpm_params(cfg, :activity)

function chemotaxis_params(cfg, fields) # chemotaxis parameters for CPM
    cpm = cfg[:rule_graph][:cpm]
    vaxes = cfg[:vaxes]
    if haskey(cpm, :chemotaxis)
        return cpm[chemotaxis]
    end

    n = types_cnt(cfg)
    vals = []

    zg = cfg[:rule_graph][:zg]
    for v in zg |> keys
        idx = findfirst(x -> x[2] == :move, zg[v])
        i = CELL_TYPES_OFFSET + vaxes[v][:id]
        if idx != nothing
            id = vaxes[zg[v][idx][1]][:id]
            weights = zeros(Float64, n)
            weights[i] = zg[v][idx][3]
            push!(vals, (weights, fields[id]))
        end
    end

    vals
end

function barrier_params(cfg) # barrier parameters
    n = types_cnt(cfg)

    barrier = fill(false, n)

    if haskey(cfg, :barrier)
        vaxes = cfg[:vaxes]
        id_barrier = CELL_TYPES_OFFSET + vaxes[cfg[:barrier][:vax]][:id]
        barrier[id_barrier] = true
    end

    barrier
end

function connectivity_params(cfg) # connectivity = constraint, that does not allow components to split
    n = types_cnt(cfg)

    connectivity = fill(false, n)

    if haskey(cfg, :connectivity)
        vaxes = cfg[:vaxes]
        vs = cfg[:connectivity][:vaxes] == :all ? vaxes |> keys : cfg[:connectivity][:vaxes]
        for v in vs
            id_connectivity = CELL_TYPES_OFFSET + vaxes[v][:id]
            connectivity[id_connectivity] = true
            end
    end

    connectivity
end

function volume_range_params(cfg) # constraint, that limits the volume of components
    n = types_cnt(cfg)

    vrange_min = fill(0, n)
    vrange_max = fill(0, n)
    has_vrange = fill(false, n)

    if haskey(cfg, :volume_range)
        vaxes = cfg[:vaxes]

        default_min = get(cfg[:volume_range], :default_min, 0)
        default_max = get(cfg[:volume_range], :default_max, cfg[:size][1]*cfg[:size][2])

        vs = (cfg[:volume_range][:vaxes] == :all
              ? map(v -> Dict(:name => v), vaxes |> keys |> collect)
              : cfg[:volume_range][:vaxes])
        for v in vs
            idx = CELL_TYPES_OFFSET + vaxes[v[:name]][:id]
            vrange_min[idx] = get(v, :min, default_min)
            vrange_max[idx] = get(v, :max, default_max)
            has_vrange[idx] = true
        end
    end

    (vrange_min, vrange_max, has_vrange)
end

function outer_action_params(cfg) # outer action (killing/adding cells)
    if !haskey(cfg, :outer_action)
        return [], nothing
    end

    regions = cfg[:outer_action][:regions]
    action = cfg[:outer_action][:action]

    regions, action
end

end # module

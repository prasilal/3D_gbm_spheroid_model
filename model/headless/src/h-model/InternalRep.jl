module InternalRep

const RREACT = 1
const RPROD = 2
const RK = 3
const RW = 4
const RR_ABSORB = 5

function nof_vaxes(cfg)
    length(keys(cfg[:vaxes]))
end

function types_cnt(cfg)
    CELL_TYPES_OFFSET + nof_vaxes(cfg)
end

function prep_vaxes_pools(cfg)
    map(cfg[:vaxes] |> collect) do (k,v)
        grid = Grid(zeros(cfg[:size]...), cfg[:torus])
        if haskey(v, :init)
            if typeof(v[:init]) <: Array
                for arr in v[:init]
                    if length(cfg[:size]) == 2
                        grid[arr[1], arr[2]] .= arr[3]
                    else
                        grid[arr[1], arr[2], arr[3]] .= arr[4]
                    end
                end
            else
                if length(cfg[:size]) == 2
                    grid[v[:init][1], v[:init][2]] .= v[:init][3]
                else
                    grid[v[:init][1], v[:init][2], v[:init][3]] .= v[:init][4]
                end
            end
        end
        grid
    end
end

function index_vaxes!(vaxes)
    for (i, k) in enumerate(keys(vaxes))
        vaxes[k][:id] = i
    end
    vaxes
end

function vax_map(vaxes)
    Dict(k => vaxes[k][:id] for k in keys(vaxes))
end

function reactions_matrix(reactions, vaxmap, col_keys = [:k, :w, :r_absorb])
    map(reactions) do r
        vcat(
            [map(v -> (v[1], vaxmap[v[2]]), r[:react]) |> collect,
             map(v -> (v[1], vaxmap[v[2]]), r[:prod]) |> collect],
            map(k -> r[k], col_keys)
        )
    end
end

function conver_zg_edge_to_idxs(edge, vaxmap)
    if edge[1] == :divide
        return tuple(edge[1:2]...,
                     (map(edge[3:end]) do r
                          (vaxmap[r[1]], r[2])
                      end)...)
    end
    edge
end

function convert_zg_to_idxs(zg, vaxmap)
    Dict([vaxmap[k],
          map(zg[k]) do es
              tuple(vaxmap[es[1]], conver_zg_edge_to_idxs(es[2:end],vaxmap)...)
          end] for k in keys(zg))
end

function convert_state_to_idxs(state, vaxmap)
    Dict((kv[1] == :cum_state ? :cum_state => vaxmap[kv[2]] : kv) for kv in state)
end

function convert_idx_to_state(state_idx, sim)
    vaxmap = sim.irep[:vaxmap]
    for (vax_name, vax_idx) in vaxmap
        if vax_idx == state_idx
            return vax_name
        end
    end
end

cum_state_to_idx(cum_state) = cum_state + CELL_TYPES_OFFSET


end #module
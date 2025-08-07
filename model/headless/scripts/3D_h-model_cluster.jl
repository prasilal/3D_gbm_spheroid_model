using Colors
const RGBAf0 = RGBA{Float32}
using MLStyle
using Random
using Distributions
using Images
using GeometryBasics
using Meshing
using JLD2


include("../src/rf.jl")
include("../src/rf_cpm.jl")
include("../src/rf_cpm/Vizs.jl")
include("../src/rf_plot.jl")

## Model

const BG_TYPE = BG_KIND # background type (1)
const CELL_TYPE_DEAD = BG_TYPE + 1 # dead cell type
const CELL_TYPES_START = CELL_TYPE_DEAD + 1 # start of cell types
const CELL_TYPES_OFFSET = CELL_TYPE_DEAD # offset for cell types

mutable struct Cell
    id :: Int
    receptors :: Vector{Float64} # receptors for each vax
    state :: Dict{Symbol, Any}
end

mutable struct Sim{N}
    cells :: Vector{Cell}
    cell_map :: Dict{Int, Cell} # map of cell id to cell
    vaxes :: Vector{Grid{Float64, N}}
    model :: GridModel 
    cfg :: Dict # configuration
    irep :: Dict{Symbol, Any} # internal representation
    in_vaxes :: Vector{Grid{Float64, N}}
end

## Conversion from CFG to internal representation

const RREACT = 1 # reactants
const RPROD = 2 # products
const RK = 3 # reaction rate
const RW = 4 # reaction weight
const RR_ABSORB = 5 # reaction absorb

function nof_vaxes(cfg) # number of vaxes
    cfg[:vaxes] |> keys |> collect |> length
end

function types_cnt(cfg) # number of types of cells
    CELL_TYPES_OFFSET + nof_vaxes(cfg)
end

function prep_vaxes_pools(cfg) 
    map(cfg[:vaxes] |> collect) do (k,v)
        grid = Grid(zeros(cfg[:size]...), cfg[:torus])
        if haskey(v, :init)
            if typeof(v[:init]) <: Array # multiple initial values
                for arr in v[:init]
                    if (length(cfg[:size]) == 2)
                        grid[arr[1],arr[2]] .= arr[3]
                    else
                        grid[arr[1],arr[2],arr[3]] .= arr[4]
                    end
                end
            else
                if (length(cfg[:size]) == 2)
                    grid[v[:init][1],v[:init][2]] .= v[:init][3] # single initial value
                else
                    grid[v[:init][1],v[:init][2], v[:init][3]] .= v[:init][4]
                end
            end
        end
        grid
    end
end

function index_vaxes!(vaxes) 
    for (i, k) in vaxes |> keys |> enumerate
        vaxes[k][:id] = i
    end

    vaxes
end

function vax_map(vaxes) # map vax names to their indeces
    Dict(k => vaxes[k][:id] for k in keys(vaxes))
end

function reactions_matrix(reactions, vaxmap, col_keys = [:k, :w, :r_absorb]) # convert reactions to matrix
    map(reactions) do r
        vcat([map(v -> (v[1], vaxmap[v[2]]), r[:react]) |> collect,
              map(v -> (v[1], vaxmap[v[2]]), r[:prod]) |> collect,],
             map(k -> r[k], col_keys))
    end
end

function conver_zg_edge_to_idxs(edge, vaxmap) # convert edge to indexes
    if edge[1] == :divide
        return tuple(edge[1:2]...,
                     (map(edge[3:end]) do r
                             (vaxmap[r[1]], r[2])
                         end)...)
    end

    edge
end

function convert_zg_to_idxs(zg, vaxmap) # convert zg to indeces
    Dict([vaxmap[k],
          map(zg[k]) do es
              tuple(vaxmap[es[1]], conver_zg_edge_to_idxs(es[2:end],vaxmap)...)
          end] for k in keys(zg))
end

function convert_state_to_idxs(state, vaxmap) # convert state to indexes
    Dict((kv[1] == :cum_state ? :cum_state => vaxmap[kv[2]] : kv) for kv in state)
end

function convert_idx_to_state(state_idx, sim) # convert index to state
    vaxmap = sim.irep[:vaxmap]
    for (vax_name, vax_idx) in vaxmap
        if vax_idx == state_idx
            return vax_name
        end
    end
end

## Preprocessing of CPM configuration 

function adhesion_params(cfg) # adhesion parameters
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

function cpm_params(cfg, key) # creates parameters for CPM (weights, values)
    cpm = cfg[:rule_graph][:cpm]
    vaxes = cfg[:vaxes]
    if haskey(cpm, key)
        return cpm[key]
    end

    n = types_cnt(cfg)
    w = zeros(Float64, n) # weights
    val = zeros(Float64, n) # values

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


## Simulation constructor

function mk_sim!(cfg, rules = Nothing)
    index_vaxes!(cfg[:vaxes])
    vaxes = prep_vaxes_pools(cfg)

    vaxmap = vax_map(cfg[:vaxes])

    internal = Dict(
        :vaxmap => vaxmap,
        :vaxes => sort(values(cfg[:vaxes]) |> collect, by = v -> v[:id]),
        :reactions => reactions_matrix(cfg[:reactions], vaxmap),
        :zg => convert_zg_to_idxs(cfg[:rule_graph][:zg], vaxmap)
    )

    params = params_cfg(
        T = get(cfg[:rule_graph][:cpm], :T, 20),
        adhesion = AdhesionParams(adhesion_params(cfg)),
        volume = VolumeParams(volume_params(cfg)...),
        perimeter = PerimeterParams(perimeter_params(cfg)...),
        chemotaxis = map(p -> ChemotaxisParams(p...), chemotaxis_params(cfg, vaxes)),
        activity = ActivityParams(Val{:geom}(), activity_params(cfg)...),
        barrier = BarrierParams(barrier_params(cfg)),
        connectivity = ConnectivityParams(connectivity_params(cfg)),
        hard_volume_range = HardVolumeRangeParams(volume_range_params(cfg)...),
        outer_action = OuterActionParams(outer_action_params(cfg)...)
    )

    rules_cfg = mk_cfg(post_mcs_listeners = [partial(rules, cfg)])
    full_cfg = merge(StatsCfg, VolumeCfg, AdhesionCfg, PerimeterCfg, ChemotaxisCfg, ActivityCfg, BarrierCfg,
                     LocalConnectivityCfg, HardVolumeRangeCfg, OuterActionCfg,
                     params, rules_cfg)
    cpm = make_preinit_cpm(cfg[:size], cfg[:seed], cfg = full_cfg, is_torus = all(identity, cfg[:torus]))

    in_vaxes = [Grid(copy(v.x), v.torus) for v in vaxes]
    cfg[:sim] = Sim(Cell[], Dict{Int, Cell}(), vaxes, cpm, cfg, internal, in_vaxes)
    cfg[:sim]
end

function mk_barrier!(cfg) # create barrier
    if !haskey(cfg, :barrier)
        return
    end

    vaxes = cfg[:vaxes]
    model = cfg[:sim].model
    id_barrier = CELL_TYPES_OFFSET + vaxes[cfg[:barrier][:vax]][:id]
    id = make_new_cell_id!(model, id_barrier)
    cfg[:barrier][:id] = id
    # cfg.sim.model.grid[cfg[:barrier][:pos]] .= id
    for idx in cfg[:barrier][:pos]
        setpix!(model, idx, id)
    end
end

function mk_cell!(sim, state, pos, init_receptors)
    id = seed_cell_at!(sim.model, cum_state_to_idx(state[:cum_state]), pos)
    cell = Cell(id, init_receptors, state)
    push!(sim.cells, cell)
    sim.cell_map[id] = cell
end

function init_receptors(receptors, vaxes)
    rs = zeros(Float64, vaxes |> keys |> length)
    for (r, v) in receptors
        rs[vaxes[r][:id]] = v
    end
    rs
end

function mk_cell_from_cfg!(sim, vaxes, c)
    state = convert_state_to_idxs(c[:state], sim.irep[:vaxmap])
    mk_cell!(sim, state, CartesianIndex(c[:init_pos]), init_receptors(c[:receptors], vaxes))
end

function mk_cells!(sim)
    cfg = sim.cfg

    vaxes = cfg[:vaxes]
    for c in cfg[:cells]
        mk_cell_from_cfg!(sim, vaxes, c)
        # state = convert_state_to_idxs(c[:state], sim.irep[:vaxmap])
        # mk_cell!(sim, state, CartesianIndex(c[:init_pos]), init_receptors(c[:receptors], vaxes))
    end
end

## Vax processes

function perform_reaction!(out_vaxes, in_vaxes, reaction, cfg)
    k = reaction[RK]
    if k == 0 return end

    r = ones(cfg[:size])

    for (n, v) in reaction[RREACT]
        @. r *= in_vaxes[v] ^ n
    end

    for (n, v) in reaction[RPROD]
        @. out_vaxes[v] += n * r * k
    end

    for (n, v) in reaction[RREACT]
        @. out_vaxes[v] -= n * r * k
    end
end

const LAPLACIAN_KERNEL_3D = OffsetArray(masked_laplacian(3, complex = true), -1:1, -1:1, -1:1)

function perform_diffusion!(sim, cfg, idx, D, barrier_idxs, dst, M, pixels_by_cell)
    if D == 0 return end

    field = sim.vaxes[idx]
    vaxes_cfg = cfg[:vaxes]
    field_id = idx

    known_ids = sim.cell_map  

    fill!(M, D)

    for (cid, pix_list) in pixels_by_cell
        if !haskey(known_ids, cid)
            continue
        end

        cell = known_ids[cid]
        cum_state = convert_idx_to_state(cell.state[:cum_state], sim)

        raw_entry = get(get(vaxes_cfg, cum_state, Dict()), :diffusion_multipliers, 1.0)
        multiplier = isa(raw_entry, Dict) ? get(raw_entry, field_id, 1.0) : raw_entry

        for idx in pix_list
            M[idx] = D * multiplier
        end
    end

    if barrier_idxs !== nothing
        M[barrier_idxs] .= 0
    end

    fill!(dst, 0.0)
    masked_manifold_smoothen!(dst, field.x, M, LAPLACIAN_KERNEL_3D)
    field.x .+= dst
end

function perform_decay!(field, d)
    field.x .*= 1.0 - d
end

function correct_vals!(field) # correct values to be non-negative
    imin = first(CartesianIndices(field.x))
    imax = last(CartesianIndices(field.x))

    @inbounds @simd for i in imin:imax
        if field.x[i] < 0
            field.x[i] = 0
        end
    end
end


function absorb_vaxes!(vax, out_vaxes, receptors, reaction, fields, cells_border_pixels) # absorbtion of vaxes by cells
    lp = 1.0
    for (n, v) in reaction[RREACT]
        lp *= (v == vax ? receptors[v] : out_vaxes[v]) ^ n
    end

    lp *= reaction[RK]
    l = length(cells_border_pixels)

    for (n, v) in reaction[RREACT]
        if v == vax
            continue
        end

        fields[v].x[cells_border_pixels] .-= n*lp/l
    end
end

function produce_vax!(field, cellpixels, val) 
    field.x[cellpixels] .+= val / length(cellpixels)
end

function produce_r_vax!(cell, edge) # produce receptors by cells
    idx = edge[1]
    cell.receptors[idx] += (get(edge, 4, 0) > 0
                            ? edge[3] * get(edge, 4, 0) * cell.state[:cum_state_weight]
                            : edge[3])
    if cell.receptors[idx] < 0
        cell.receptors[idx] = 0
    end
end

function produce_v_vax!(cell, edge, sim, fields) #produce vax by cells
    vaxes = sim.cfg[:vaxes]
    cellpixels = get(get_stat(sim.model, :border_pixels_by_cell), cell.id, [])

    field = fields[edge[1]]
    vol = (get(edge, 4, 0) > 0
           ? get(edge, 4, 0) * cell.state[:cum_state_weight] * edge[3]
           : edge[3])
    produce_vax!(field, cellpixels, vol)
    for idx in cellpixels
        if field.x[idx] < 0
            field.x[idx] = 0
        end
    end
end

## Vax functionality

function reactions_for(reactions, vax) # finds all reactions for vax
    filter(reactions) do r
        findfirst(x -> x[2] == vax, r[RREACT]) != nothing
    end
end

function get_outer_vaxes(vaxes, cellpixels) # get concentration of vaxes outside of cells
    map(vaxes) do vax
        vax[cellpixels] |> sum
    end
end

function get_bound_vaxes(nof_vaxes, sim, cell_border_list, cells_border_pixels) # what receptors are on neighboring cells
    bound_vaxes = zeros(Float64, nof_vaxes)
    for id in (@>> cell_border_list keys collect filter(isneqv(BGID)))
        if haskey(sim.cell_map, id)
            cell = sim.cell_map[id]
            border_weight = cell_border_list[id] / length(cells_border_pixels[id])
            bound_vaxes .+= cell.receptors .* border_weight
        end
    end
    bound_vaxes
end

## Cell functionality

cum_state_to_idx(cum_state) = cum_state + CELL_TYPES_OFFSET # convert cum state to index

function cur_cum_state(sim, the_cell, in_vaxes) # infers what the state should be 
    cfg = sim.cfg
    reactions = sim.irep[:reactions]
    id = the_cell.id

    cells_border_pixels = get_stat(sim.model, :border_pixels_by_cell)
    cell_border_list = get(get_stat(sim.model, :cell_neighbor_list), id, Dict())

    n = cfg[:vaxes] |> length
    out_vaxes = get_outer_vaxes(in_vaxes, get(cells_border_pixels, id, []))
    bound_vaxes = get_bound_vaxes(n, sim, cell_border_list, cells_border_pixels)
    out_bound_vaxes = out_vaxes .+ bound_vaxes

    weights = zeros(Float64, n)
    for vax in 1:n
        rs = reactions_for(reactions, vax)
        wp = 0.0
        for r in rs
            if r[RW] == 0
                continue
            end

            lp = 1.0
            for (n, v) in r[RREACT]
                lp *= (v == vax ? the_cell.receptors[v] : out_bound_vaxes[v]) ^ n
            end
            wp += lp * r[RW]

            if r[RR_ABSORB]
                absorb_vaxes!(vax, out_vaxes, the_cell.receptors,
                              r, sim.vaxes, get(cells_border_pixels, id, []))
            end
        end
        weights[vax] = wp
    end

    idx = argmax(weights)

    if iszero(weights) || weights[idx] < cfg[:rule_graph][:min_weight] || weights[idx] == 0.0
        the_cell.state[:cum_state], the_cell.state[:cum_state_weight]
    else
        idx, weights[idx]
    end
end

function cell_divide!(cell, edge, sim) # cell division
    rs = cell.receptors |> sum
    max_rs = typeof(edge[3]) <: Number ? edge[3] : edge[3][1]
    max_rs_div = typeof(edge[3]) <: Number ? 0.5 : edge[3][2]

    if rs < max_rs
        return nothing
    end

    println("diving cell $(cell.id). num of receptors $rs threshold $(edge[3])")

    nid = divide_cell!(sim.model, cell.id)
    new_cell = Cell(nid, deepcopy(cell.receptors), deepcopy(cell.state))
    mult = fill(0.5, length(cell.receptors))

    for (r, m) in edge[4:end]
        mult[r] = m
    end

    cell.receptors .*= mult

    if sum(cell.receptors) >= max_rs
        cell.receptors .*= max_rs_div
    end

    @. mult = 1.0 - mult
    new_cell.receptors .*= mult

    if sum(new_cell.receptors) >= max_rs
        new_cell.receptors .*= max_rs_div
    end

    new_cell.state[:cum_state] = edge[1]

    new_cell
end

function necrosis!(cell, sim, vaxes) # spills vaxes from dead cells
    cellpixels = get(get_stat(sim.model, :pixels_by_cell), cell.id, [])

    for idx in 1:length(cell.receptors)
        produce_vax!(vaxes[idx], cellpixels, cell.receptors[idx])
    end
end

function hard_kill!(cell, sim) # deletes cell (apoptosis)
    kill_cell!(sim.model, cell.id)
end

function hard_kill_n_del!(cid, sim)
    cell = sim.cell_map[cid]
    hard_kill!(cell, sim)
    delete!(sim.cell_map, cell.id)
    deleteat!(sim.cells, findfirst(c -> c.id == cell.id, sim.cells))
end


function soft_kill!(cell, sim) # cell still occupies space, but is dead
    set_cell_kind!(sim.model, cell.id, CELL_TYPE_DEAD)
end

function cell_kill!(cell, edge, sim, vaxes) 
    if rand() > edge[3]
        return false
    end

    println("killing cell $(cell.id) by $edge")

    @match edge[4] begin
        :apoptosis => nothing
        :necrosis => necrosis!(cell, sim, vaxes)
    end

    @match edge[5] begin
        :hard => hard_kill!(cell, sim)
        :soft => soft_kill!(cell, sim)
    end

    true
end

function update_by_zg!(sim, cell, vaxes) # updates cell by zg
    killed, new_cell = false, nothing

    for edge in sim.irep[:zg][cell.state[:cum_state]]
        action = edge[2]
        if action == :prod_r
            produce_r_vax!(cell, edge)
        elseif action == :prod_v
            produce_v_vax!(cell, edge, sim, vaxes)
        elseif action == :divide
            cell_divide_fn = get(get(sim.cfg, :rule_graph, Dict()), :division_hook, cell_divide!)
            new_cell = cell_divide_fn(cell, edge, sim)
        elseif action == :kill
            killed = cell_kill!(cell, edge, sim, vaxes)
        end
    end

    res = killed ? [] : [cell]
    if new_cell != nothing
        push!(res, new_cell)
        sim.cell_map[new_cell.id] = new_cell
    end

    res
end

function get_resting_time(graph, zg, cstate) # time until next state change
    zg_actions = zg[cstate]
    rt_pos = findfirst(x -> x[2] == :resting_time, zg_actions)

    if rt_pos == nothing
        return graph[:resting_time]
    end

    zg_actions[rt_pos][3]
end

function update_cum_state(sim, cell, in_vaxes) # updates cell state
    if cell.state[:resting_time] > 0
        cell.state[:resting_time] -= 1
        return
    end

    cstate, cweight = cur_cum_state(sim, cell, in_vaxes)
    if cell.state[:cum_state] != cstate
        cell.state[:resting_time] = get_resting_time(sim.cfg[:rule_graph], sim.irep[:zg], cstate)

        cell_state_update = get(sim.cfg, :cell_state_update, nothing)
        if cell_state_update != nothing
            for f in get(cell_state_update, :hooks, [])
                f(sim, cell, cstate, cweight, get(cell_state_update, :state, nothing))
            end
        end

        if get(get(sim.cfg[:runtime], :debug, Dict()), :show_cell_state_change, true) == true
            println("Update cum state for cell $(cell.id): old($(convert_idx_to_state(cell.state[:cum_state], sim)), $(cell.state[:cum_state_weight])); new($(convert_idx_to_state(cstate, sim)), $cweight)")
        end
    end
    cell.state[:cum_state] = cstate
    cell.state[:cum_state_weight] = cweight
    set_cell_kind!(sim.model, cell.id, cum_state_to_idx(cstate))
end

function update_r_vaxes(sim, cell) # updates receptors (decay)
    vaxes = sim.irep[:vaxes]

    for (idx, vinfo) in enumerate(vaxes)
        cell.receptors[idx] *= 1 - vinfo[:rd]
    end
end

#### Rules
function absorb_vaxes_cell!(sim::Sim{3}, cell::Cell, in_vaxes::Vector{Grid{Float64, 3}})
    return  # do nothing
end


function rules(cfg, model :: GridModel{T, N, CPMState{Tdc, Trnd, N}, CPMCfg}) where {T, Tdc, N, Trnd}
    if !get(cfg, :rules_enabled, true) return end

    sim = cfg[:sim]
    in_vaxes = sim.in_vaxes

    for j in 1:length(sim.vaxes)
        copyto!(in_vaxes[j].x, sim.vaxes[j].x)
    end

    barrier_idxs = haskey(cfg, :barrier) ? cfg[:barrier][:pos] : nothing
    update_cum_state_fn = get(get(sim.cfg, :cell_state_update, Dict()), :transition_fn, update_cum_state)
    absorb_vaxes_fn = get(get(sim.cfg, :cell_state_update, Dict()), :absorb_vaxes_for_cell, absorb_vaxes_cell!)

    
    Threads.@threads for i in eachindex(sim.cells)
        update_cum_state_fn(sim, sim.cells[i], in_vaxes)
    end



    cells = Vector{Cell}()
    for c in sim.cells
	absorb_vaxes_fn(sim, c, sim.vaxes)
        new_cells = update_by_zg!(sim, c, sim.vaxes)
        update_r_vaxes(sim, c)
        append!(cells, new_cells)
    end
    sim.cells = cells

    reset_stat_state!(sim.model)

    nthreads = Threads.nthreads()
    sz = size(sim.vaxes[1])
    dsts = [zeros(Float64, sz) for _ in 1:nthreads]
    Ms = [fill(1.0, sz) for _ in 1:nthreads]

    for i in 1:get(cfg, :vaxes_rule_steps, 1)
        if i > 1
            for j in 1:length(sim.vaxes)
                copyto!(in_vaxes[j].x, sim.vaxes[j].x)
            end
	end

        for r in sim.irep[:reactions]
            perform_reaction!(sim.vaxes, in_vaxes, r, cfg)
        end

        vaxes = sim.irep[:vaxes]
	pixels_by_cell = get_stat(sim.model, :pixels_by_cell)

        Threads.@threads for idx in 1:length(vaxes)
            tid = Threads.threadid()
            vinfo = vaxes[idx]
            field = sim.vaxes[idx]
            dst = dsts[tid]
            M = Ms[tid]

            perform_diffusion!(sim, cfg, idx, vinfo[:D], barrier_idxs, dst, M, pixels_by_cell)
            perform_decay!(field, vinfo[:d])
            correct_vals!(field)
        end
    end
end

#### Visualization

rc(hex) = ((hex & 0xff0000) >> 16) / 255
gc(hex) = ((hex & 0xff00) >> 8) / 255
bc(hex) = (hex & 0xff) / 255
hex2rgba(hex) = RGBA(rc(hex), gc(hex), bc(hex), 1.)

default_colors = Dict( # default colors for cells
    0 => hex2rgba(0x370617), # cell border
    1 => RGBA(1,1,1,1), # bg
    2 => hex2rgba(0x9fa6a6),
    3 => hex2rgba(0xf1f3c2),
    4 => hex2rgba(0xf3e700),
    5 => hex2rgba(0xf3d368),
    6 => hex2rgba(0xf3ba9e),
    7 => hex2rgba(0xf37b35),
    8 => hex2rgba(0xf31000),
    9 => hex2rgba(0xf3b4f1),
    10 => hex2rgba(0xf300ae),
    11 => hex2rgba(0x9a00f3),
    12 => hex2rgba(0x6100f3),
    13 => hex2rgba(0xb48bf3),
    14 => hex2rgba(0x648ff3),
    15 => hex2rgba(0xadd6f3),
    16 => hex2rgba(0x1119f3),
    17 => hex2rgba(0x1df3e8),
    18 => hex2rgba(0x20f386),
    19 => hex2rgba(0xbff3a3),
    20 => hex2rgba(0x00a876),
    21 => hex2rgba(0x961200),
    22 => hex2rgba(0xffeedb)
)

default_colors_field = Dict( # default colors for vaxes
    1 => hex2rgba(0xf1f3c2),
    2 => hex2rgba(0xf3e700),
    3 => hex2rgba(0xf3d368),
    4 => hex2rgba(0xf3ba9e),
    5 => hex2rgba(0xf37b35),
    6 => hex2rgba(0xf31000),
    7 => hex2rgba(0xf3b4f1),
    8 => hex2rgba(0xf300ae),
    9 => hex2rgba(0x9a00f3),
    10 => hex2rgba(0x6100f3),
    11 => hex2rgba(0xb48bf3),
    12 => hex2rgba(0x648ff3),
    13 => hex2rgba(0xadd6f3),
    14 => hex2rgba(0x1119f3),
    15 => hex2rgba(0x1df3e8),
    16 => hex2rgba(0x20f386),
    17 => hex2rgba(0xbff3a3),
    18 => hex2rgba(0x00a876),
    19 => hex2rgba(0x961200),
    20 => hex2rgba(0xffeedb)
)

function draw_cpm(sim;
                  cell_kind_color = default_colors,
                  num_cell_kinds = nothing,
                  show_activity = false,
                  show_fields = false,
                  field_color = default_colors_field,
                  vox = nothing)
    cpm = sim.model
    if show_fields
        for v in values(sim.cfg[:vaxes])
            if get(v,:show,false)
                # vox = draw_field_contour(sim.vaxes[v[:id]],
                #                          haskey(v,:color) ? v[:color] : field_color[v[:id]],
                #                          vox = vox, nsteps = 30)
                vox = draw_field(sim.vaxes[v[:id]],
                                 haskey(v,:color) ? v[:color] : field_color[v[:id]],
                                 vox = vox)
            end
        end
    end

    vox = draw_cells(cpm, cid -> cell_kind_color[cell_kind(cpm, cid)], vox = vox)

    num_cell_kinds = num_cell_kinds == nothing ? maximum(cpm.state.cell_to_kind) : num_cell_kinds

    for k in 2:num_cell_kinds
        draw_cell_borders(cpm, k, cell_kind_color[0], vox = vox)
        if show_activity
            draw_activity_values(cpm, k, vox = vox)
        end
    end

    vox
end

### 3D visualization helper functions - not used yet

function draw_cpm_3d(sim;
                    cell_kind_color = default_colors,
                    num_cell_kinds = nothing,
                    show_activity = false,
                    show_fields = false,
                    field_color = default_colors_field,
                    vox = nothing)
    
    # Delete old meshes
    for p in cell_meshes_obs[]
        delete!(ax1.scene, p)
    end
    cell_meshes_obs[] = Plot[]

    update_cell_colors!()

    # Draw 3D cubes for cell voxels
    cell_voxels = get_stat(sim.model, :pixels_by_cell)
    for (cid, voxels) in cell_voxels
        kind = cell_kind(sim.model, cid)
        color = get(cell_color_map, cid, RGBAf0(1, 1, 1, 1))
        for v in voxels
            x, y, z = Tuple(v)
            p = mesh!(ax1, Rect3f(Vec3f(x, y, z), Vec3f(1, 1, 1)), color=color, transparency=true)
            push!(cell_meshes_obs[], p)
        end
    end

    # Draw borders in black for each kind (from 2 onward)
    for k in 2:(num_cell_kinds === nothing ? maximum(sim.model.state.cell_to_kind) : num_cell_kinds)
        border_vox = draw_cell_borders(sim.model, k, RGBAf0(1, 1, 1, 1))  # returns RGBA array
        border_indices = findall(x -> x.alpha > 0.0, border_vox)
        for idx in border_indices
            x, y, z = Tuple(idx)
            b = mesh!(ax1, Rect3f(Vec3f(x, y, z), Vec3f(1, 1, 1)), color=RGBAf0(0, 0, 0, 1.0))
            push!(cell_meshes_obs[], b)
        end
    end


    # Update fields
    if show_fields
        for (key, volplot) in volumes
            field = sim.vaxes[sim.cfg[:vaxes][key][:id]]
            vol_data = convert(Array{Float32, 3}, field.x)
            volplot[4][] = vol_data  
        end
    end

    step_label.text[] = "Step: $(Int(round(t[])))"
    count_label.text[] = "Cell count: $(length(sim.cells))"
    
    pixels_by_cell = get_stat(sim.model, :pixels_by_cell)
    cell_volumes = [length(get(pixels_by_cell, c.id, [])) for c in sim.cells]
    
    vol_label.text[] = "Avg vol: $(round(mean(cell_volumes), digits=1))"
    max_vol_label.text[] = "Max vol: $(maximum(cell_volumes))"            
    
    receptor_vals = [sum(c.receptors) for c in sim.cells]
    receptor_label.text[] = "Avg receptor: $(round(mean(receptor_vals), digits=3))"
end

function draw_cells_3D(sim; ax=nothing)
    scene = ax === nothing ? LScene(size=(800,800)) : ax

    model = sim.model
    state = model.state
    pxs = state.extra[:pixels_by_cell]
    colors = sim.cfg[:runtime][:cell_kind_color]

    for (cid, coords) in pxs
        clr = get(colors, sim.model.state.kinds[cid], RGBA(1, 0, 0, 1.0))
        for pos in coords
            x, y, z = Tuple(pos)
            mesh!(scene, Rect3f(Vec{3, Float32}(x, y, z), Vec{3, Float32}(1, 1, 1)), color=clr, transparency=true)
        end
    end

    scene
end






### Simulation control

always_true(_) = true

function time_macro_step!(sim)
    rules_enabled = get(sim.cfg,:rules_enabled, true)
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
        save_sims_state(sim, i)
        if i % 5 == 0
            GC.gc() 
        end
        if !should_we_continue(sim)
            break
        end
    end
end

### Visualization prep

function draw_smooth_cell_mesh!(ax, pixels::Vector{CartesianIndex{3}}, grid_size::Tuple{Int, Int, Int}, color::RGBA)
    vol = zeros(Float32, grid_size)
    for p in pixels
        vol[p] = 1.0
    end

    xs = 1:grid_size[1]
    ys = 1:grid_size[2]
    zs = 1:grid_size[3]

    points, faces = isosurface(vol, MarchingCubes(iso=0.5), xs, ys, zs)

    gb_points = GeometryBasics.Point{3, Float32}.(points)
    gb_faces = GeometryBasics.TriangleFace.(faces)
    gb_mesh = GeometryBasics.Mesh(gb_points, gb_faces)

    return mesh!(ax, gb_mesh, color=color, transparency=true)
end

function get_cell_color(cell::Cell, sim::Sim)
    state_sym = convert_idx_to_state(cell.state[:cum_state], sim)
    parts = split(String(state_sym), "_")
    base = Symbol(parts[1])           
    subtype = Symbol(parts[2])        

    base_colors = Dict(
        :gbm => RGB(1.0, 0.0, 0.0),
        :msc => RGB(0.0, 0.0, 1.0),
        :gsc => RGB(0.0, 1.0, 0.0),
    )

    state_intensity = Dict(
        :p => 1.0,   
        :q => 0.6,
        :n => 0.3
    )

    base_rgb = get(base_colors, base, RGB(0.5, 0.5, 0.5))
    intensity = get(state_intensity, subtype, 0.5)

    scaled_rgb = RGB(
        base_rgb.r * intensity + (1.0 - intensity),
        base_rgb.g * intensity + (1.0 - intensity),
        base_rgb.b * intensity + (1.0 - intensity),
    )

    return RGBAf0(scaled_rgb.r, scaled_rgb.g, scaled_rgb.b, 0.9f0)
end


function prep_vis_sim(sim;
    burnin_steps = 20,
    show_3d = true,
    num_cell_kinds = nothing,
    show_activity = false,
    show_fields = false,
    cell_kind_color = default_colors,
    field_color = default_colors_field)
    return nothing, nothing
end
    






function anim_sim(sim, t; num_of_steps = 1000, fps = 1.0/30., should_we_continue = always_true)
    for i in 1:num_of_steps
        t[] += 1
        if !should_we_continue(sim)
            break;
        end
        sleep(fps)
    end
end


function record_sim(sim, t, img; filename = "simulation.mp4", num_of_steps = 100, fps = 1.0 / 30.0, should_we_continue = always_true)
    timestamps = range(0, step = fps, length = num_of_steps)
    record(img, filename, timestamps; framerate = 1.0 / fps) do _
        t[] += 1
        if !should_we_continue(sim)
            stop()
        end
    end
end


const DEFAULT_RUNTIME_PARAMS = @d(
    :burnin_steps => 20,
    :micro_steps => 1,
    :show_fields => true,
    :show_activity => true,
    :show_sim => true,
    :cell_kind_color => default_colors,
    :field_color => default_colors_field
)

function get_black_pixel_coordinates_scaled(image_path, target_width, target_height; threshold = 0.5)
    img_org = load(image_path)
    img = Gray.(img_org)

    height, width = size(img)

    black_pixel_coordinates = Set{Tuple{Int, Int}}()

    x_scale = target_width / width
    y_scale = target_height / height

    for y in 1:height
        for x in 1:width

            pixel_value = img[y, x]

            if pixel_value <= Gray(threshold)

                scaled_x = round(Int, x * x_scale)
                scaled_y = round(Int, y * y_scale)

                push!(black_pixel_coordinates, (scaled_x+1, scaled_y+1))
            end
        end
    end

    return [CartesianIndex(idx) for idx in black_pixel_coordinates]
end

function init_sim(cfg)
    sim = mk_sim!(cfg, rules)
    mk_barrier!(cfg)
    mk_cells!(sim)

    sim.cell_map = Dict(c.id => c for c in sim.cells)

    rules_enabled = get(cfg, :rules_enabled, true)

    img, t = nothing, nothing

    cfg[:runtime] = merge(DEFAULT_RUNTIME_PARAMS, get(cfg,:runtime,@d()))

    cfg[:rules_enabled] = false
    if cfg[:runtime][:show_sim]
        img,t = prep_vis_sim(sim,
                             burnin_steps = cfg[:runtime][:burnin_steps],
                             show_3d= cfg[:runtime][:show_3d],
                             num_cell_kinds = CELL_TYPES_OFFSET+length(sim.vaxes),
                             show_fields = cfg[:runtime][:show_fields],
                             show_activity = cfg[:runtime][:show_activity],
                             cell_kind_color = cfg[:runtime][:cell_kind_color],
                             field_color = cfg[:runtime][:field_color])
    else
        burnin!(sim, cfg[:runtime][:burnin_steps])
    end
    cfg[:rules_enabled] = rules_enabled

    @d(
        :sim => sim,
        :img => img,
        :t => t
    )
end

function simulate(sim_desc; num_of_steps = 1, record = false, filename = "simulation.mp4", fps = 1.0/30., should_we_continue = always_true)
    if sim_desc[:sim].cfg[:runtime][:show_sim]
        if record
            record_sim(sim_desc[:sim], sim_desc[:t], sim_desc[:img], filename = filename, num_of_steps = num_of_steps, fps = fps, should_we_continue = should_we_continue)
        else
            display(sim_desc[:img])
            anim_sim(sim_desc[:sim], sim_desc[:t], num_of_steps = num_of_steps, fps = fps, should_we_continue = should_we_continue)
        end
    else
        burnin!(sim_desc[:sim], num_of_steps, should_we_continue = should_we_continue)
    end
    #sim_desc[:sim]
end




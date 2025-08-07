using Colors
const RGBAf0 = RGBA{Float32}
using GLMakie
using MLStyle
using Random
using Distributions
using Images
using GeometryBasics
using Meshing


include("../src/rf.jl")
include("../src/rf_cpm.jl")
include("../src/rf_cpm/Vizs.jl")
include("../src/rf_plot.jl")

## Model

const BG_TYPE = BG_KIND # background type (1)
const CELL_TYPE_DEAD = BG_TYPE + 1 # dead cell type
const CELL_TYPES_START = CELL_TYPE_DEAD + 1 # start of cell types
const CELL_TYPES_OFFSET = CELL_TYPE_DEAD # offset for cell types


"""
    Cell

A biological cell in the simulation, represented as an object with:
- `id`: Unique integer identifier, corresponding to the CPM label
- `receptors`: Vector of receptor levels (one per chemical field, or "vax")
- `state`: Dictionary of cell state attributes, including:
    - `:cum_state` — the current cell identity (e.g. `:gbm_p`)
    - `:cum_state_weight` — internal weight from field sensing
    - `:resting_time` — number of steps before a state change is allowed
    - any other metadata needed for rules

This structure holds all per-cell biological and regulatory information.
"""
mutable struct Cell
    id :: Int
    receptors :: Vector{Float64} # receptors for each vax
    state :: Dict{Symbol, Any}
end


"""
    Sim{N}

The core simulation object representing a Cellular Potts Model with chemical fields and active cells.

Fields:
- `cells`: List of all live `Cell` objects
- `cell_map`: Dict mapping cell IDs to `Cell` objects for fast lookup
- `vaxes`: List of `Grid{Float64, N}` voxel fields (chemical concentrations) for each defined variable
- `model`: The `GridModel` holding CPM lattice state, energy, and boundary information
- `cfg`: Full configuration dictionary for the simulation, containing rules, CPM parameters, field setup, etc.
- `irep`: Internal representations (e.g., index maps for reaction targets, precompiled rule graph, etc.)
- `in_vaxes`: Snapshot of previous field state used for reaction/diffusion calculations

This struct encapsulates the evolving simulation state during runtime.
"""
mutable struct Sim{N}
    cells :: Vector{Cell}
    cell_map :: Dict{Int, Cell}
    vaxes :: Vector{Grid{Float64, N}}
    model :: GridModel
    cfg :: Dict
    irep :: Dict{Symbol, Any}
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

"""
    mk_sim!(cfg::Dict, rules = Nothing) -> Sim

Constructs the full simulation object (`Sim`) from a configuration dictionary `cfg`.

This function:
- Converts the user-defined configuration into internal representations (`:vaxmap`, `:zg`, `:reactions`)
- Creates the CPM model with proper constraints (adhesion, volume, chemotaxis, etc.)
- Prepares the voxel fields for each chemical variable (vaxes)
- Initializes data structures for receptor diffusion, field interactions, and rule-based updates

If a `rules` function is provided, it is registered as a `post_mcs_listener` for the simulation steps.

Returns the `Sim` object stored in `cfg[:sim]`. This function is typically followed by `mk_cells!` and `mk_barrier!`.
"""
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


"""
    mk_barrier!(cfg::Dict)

Creates a static barrier region within the CPM lattice based on `cfg[:barrier]`.

The barrier setup includes:
- `:vax`: The symbolic name of the vax associated with the barrier (used to determine its cell kind index)
- `:pos`: A vector of `CartesianIndex` positions where the barrier should be placed

Behavior:
- A new CPM cell is created and placed at each specified voxel in `:pos`
- The cell is assigned a unique ID and given the corresponding cell kind based on the vax
- The barrier ID is stored in `cfg[:barrier][:id]` for later reference (e.g., to block diffusion)

This function must be called **after** `mk_sim!` but **before** simulation starts.

Note: Barrier voxels are not added to the list of active cells and do not participate in cell updates.
"""
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


"""
    mk_cell!(sim::Sim, state::Dict{Symbol, Any}, pos::CartesianIndex, init_receptors::Vector{Float64})

Places a single cell into the CPM lattice and registers it in the simulation.

Arguments:
- `sim`: The current `Sim` object
- `state`: A dictionary describing the cell's state (must include `:cum_state`)
- `pos`: Initial lattice position (as a `CartesianIndex`) where the cell will be seeded
- `init_receptors`: Initial receptor levels for all chemical fields (vaxes)

Behavior:
- Creates a new CPM cell at the specified position using `cum_state_to_idx` to set the cell kind
- Instantiates a `Cell` object with the given state and receptors
- Registers the cell in `sim.cell_map` and appends it to `sim.cells`

This function is typically called internally by `mk_cells!` when initializing the simulation.
"""
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


"""
    mk_cells!(sim::Sim)

Initializes and places all cells in the simulation according to the configuration `sim.cfg[:cells]`.

For each cell:
- Converts symbolic `:state` into internal indexed form using `:vaxmap`
- Initializes receptor vector based on specified values in `:receptors`
- Seeds the cell at the provided lattice position (`:init_pos`)
- Registers the cell in `sim.cell_map` and `sim.cells`

This function must be called after `mk_sim!` and before simulation steps begin.

Note: Positions must be valid within the grid and non-overlapping.
"""
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

"""
    perform_reaction!(out_vaxes, in_vaxes, reaction, cfg)

Applies a single chemical reaction to the voxel fields during a simulation step.

Arguments:
- `out_vaxes`: Output fields (modified in-place)
- `in_vaxes`: Input snapshot of field values (immutable during this step)
- `reaction`: A pre-parsed reaction tuple from `reactions_matrix`, of the form:
    `[reactants, products, k, w, absorb_flag]`
- `cfg`: Global configuration dictionary

Reaction behavior:
- Calculates a field-wide reaction term `r` based on the product of reactant concentrations at each voxel
- Updates each product field by adding `+n·r·k` and each reactant by subtracting `−n·r·k`
- No receptor or per-cell contributions are included — this is purely a spatial field-level reaction

No return value; modifies `out_vaxes` in-place.

This function is called multiple times per step inside the `rules` loop.
"""
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


"""
    perform_diffusion!(sim, cfg, idx::Int, D::Float64, barrier_idxs, dst, M, pixels_by_cell)

Performs masked anisotropic diffusion of a single voxel field with optional spatial heterogeneity.

Arguments:
- `sim`: Simulation object
- `cfg`: Configuration dictionary
- `idx`: Index of the voxel field to update
- `D`: Baseline diffusion coefficient
- `barrier_idxs`: Optional list of Cartesian indices that block diffusion (e.g., for physical barriers)
- `dst`: Temporary buffer for Laplacian output
- `M`: Per-voxel diffusion multiplier map (modified during execution)
- `pixels_by_cell`: Dict mapping cell IDs to their occupied voxels

Diffusion behavior:
- Resets multiplier `M` to `D`, then scales it per voxel depending on:
  - The cell occupying that voxel (via `pixels_by_cell`)
  - The cell’s `cum_state` and corresponding field-specific diffusion multiplier from the config
- Zeros out diffusion in barrier regions if `barrier_idxs` is provided
- Applies a masked Laplacian (`masked_manifold_smoothen!`) to `field.x` using kernel `LAPLACIAN_KERNEL_3D`

No return value; modifies the field in-place.

This method allows spatially heterogeneous and state-dependent diffusion, as seen in real tumor tissue.
"""
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


"""
    cell_divide!(cell::Cell, edge::Tuple, sim::Sim) -> Union{Cell, Nothing}

Performs mitotic division of a single cell if it meets the receptor threshold defined by `edge[3]`.

Division logic:
- Cell must have total receptor sum above a minimum threshold
- Receptors are split between parent and daughter using `edge[4:end]` as optional per-receptor multipliers
- If either cell still exceeds the max threshold post-division, receptors are scaled down again
- The daughter inherits the parent’s state, with its `:cum_state` overwritten by `edge[1]`
- Registers the daughter in `sim.cell_map`, but does not add it to `sim.cells` (caller must do that)

Returns the new `Cell` object if division occurred, or `nothing` otherwise.
"""
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

"""
    update_by_zg!(sim::Sim, cell::Cell, vaxes::Vector{Grid{Float64}}) -> Vector{Cell}

Applies rule-graph (`:zg`) actions to a single cell based on its current `:cum_state`.

For each edge in the rule graph:
- Executes actions such as:
  - `:prod_r` - receptor production
  - `:prod_v` - secretion of chemical fields (vaxes)
  - `:divide` - triggers cell division via `cell_divide!`
  - `:kill` - induces cell death via apoptosis or necrosis
- Updates cell state in-place or generates a new daughter cell if division occurs

Returns a vector of resulting live cells (`[original]` or `[original, daughter]`). If the cell was killed, returns an empty vector.

This function is typically called once per cell during each macrostep.
"""
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
function absorb_vaxes_cell!(sim::Sim{3}, cell::Cell, in_vaxes::Vector{Grid{Float64, 3}})
    return  
end


"""
    rules(cfg, model::GridModel)

Applies the full rule set for a single macrostep update in the simulation.

This includes:
- Cell state transitions based on field sensing and logic hooks
- Cell actions defined in the rule graph (`:zg`), such as receptor production, vax secretion, division, or death
- Receptor decay for each cell
- Chemical field reactions, diffusion, and decay across multiple substeps (`cfg[:vaxes_rule_steps]`)
- Optional constraints like barriers and field-specific diffusion multipliers

This function is intended to be used as the simulation "rule engine" and is passed as a callback into `mk_sim!`.
"""
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

### Simulation control

always_true(_) = true


"""
    time_macro_step!(sim::Sim)

Performs a single **macrostep** in the simulation, consisting of one full CPM update and optional rule application.

Behavior:
- Executes `cfg[:runtime][:micro_steps] - 1` raw CPM updates with `rules_enabled = false`
- Then performs one final update with `rules_enabled` restored (allowing biological rules like division, death, or secretion to apply)
- This mechanism separates physical dynamics (movement, adhesion) from biological behaviors (defined in `rules`)

Intended for use in interactive stepping or animation (`prep_vis_sim`, `anim_sim`, etc.).

Note: Assumes the `rules` function has been set up as a `post_mcs_listener` in the CPM model.
"""
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


"""
    burnin!(sim::Sim, steps::Int; should_we_continue = always_true)

Runs a series of CPM macrosteps **without visualization or recording**, typically to stabilize cell positions before measurement.

Arguments:
- `sim`: The simulation object
- `steps`: Number of macrosteps to simulate
- `should_we_continue`: Optional callback `(sim) -> Bool` to support early termination

Behavior:
- Repeatedly calls `time_macro_step!` for the specified number of steps
- Used during initialization (e.g., `prep_vis_sim`) to allow cells to settle before tracking dynamics

Returns nothing. Operates in-place.
"""
function burnin!(sim, steps; should_we_continue = always_true)
    for i in 1:steps
        time_macro_step!(sim)
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


"""
    prep_vis_sim(sim::Sim; burnin_steps::Int = 20, show_3d::Bool = true, ... ) -> (Figure, Observable)

Prepares the interactive 3D or 2D visualization for a running simulation using GLMakie.

Arguments:
- `burnin_steps`: Number of initial macrosteps before the visualization starts (e.g., to equilibrate cell positions)
- `show_3d`: If true, sets up 3D visualization with fields, stats, and controls; otherwise uses simpler 2D view
- `num_cell_kinds`: Number of distinct cell kinds (usually inferred automatically)
- `show_activity`: If true, overlays activity values
- `show_fields`: If true, enables display of chemical fields (vaxes)
- `cell_kind_color`, `field_color`: Color maps for cells and fields

Returns:
- `fig`: A `Makie.Figure` instance containing all interactive plots and UI elements
- `t`: A `Makie.Observable` time object that drives updates with each macrostep

This function is internally called by `init_sim` if `cfg[:runtime][:show_sim]` is true.
"""
function prep_vis_sim(sim;
    burnin_steps = 20,
    show_3d = true,
    num_cell_kinds = nothing,
    show_activity = false,
    show_fields = false,
    cell_kind_color = default_colors,
    field_color = default_colors_field)

    burnin!(sim, burnin_steps)
    t = Observable(0.0)

    cell_color_map = Dict{Int, RGBA{Float32}}()
    
    update_cell_colors!() = begin
        ids = [c.id for c in sim.cells]
        raw_colors = distinguishable_colors(length(ids))
        for (i, id) in enumerate(ids)
            rgb = raw_colors[i]
            cell_color_map[id] = RGBAf0(rgb.r, rgb.g, rgb.b, 0.9f0)
        end
    end
    
    update_cell_colors!() 

    if show_3d
        fig = Figure(size = (1200, 800))  # Big canvas
        ax1 = Axis3(fig[1, 1], title = "Cells", aspect = :data)
        ax2 = Axis3(fig[1, 2], title = "Vaxes", aspect = :data)
        
        # Stat table
        stats_ax = fig[2, 1:2] = GridLayout()
        
        # Stat labels
        step_label = Label(stats_ax[1, 1], "Step: 0", halign = :left)
        count_label = Label(stats_ax[1, 2], "Cell count: 0", halign = :left)
        vol_label = Label(stats_ax[1, 3], "Avg vol: 0", halign = :left)
        max_vol_label = Label(stats_ax[1, 4], "Max vol: 0", halign = :left)
        state_label = Label(stats_ax[2, 1:4], "States: ", halign = :left)

        # Showable vaxes
        shown_vaxes = show_fields ? [k for (k, v) in sim.cfg[:vaxes] if get(v, :show, false)] : Symbol[]
        selected_index = Observable(1)
        selected_vax = lift(selected_index) do i
            isempty(shown_vaxes) ? :none : shown_vaxes[i]
        end

        # Button + label for switching fields
        vax_label = Label(stats_ax[3, 1], "Field: $(selected_vax[])")
        next_button = Button(stats_ax[3, 2], label = "Next Field →")
        on(next_button.clicks) do _
            if !isempty(shown_vaxes)
                selected_index[] = selected_index[] % length(shown_vaxes) + 1
                vax_label.text[] = "Field: $(selected_vax[])"
            end
        end

        # Grid limits
        size_x, size_y, size_z = size(sim.model.grid)
        limits!(ax1, 0, size_x, 0, size_y, 0, size_z)

        # Field volume plots
        volumes = Dict{Symbol, Any}()
        if show_fields
            for key in shown_vaxes
                field = sim.vaxes[sim.cfg[:vaxes][key][:id]]
                vol = volume!(ax2, field.x; algorithm = :mip, colormap = :viridis, transparency = true)
                volumes[key] = vol
            end
        end

        cell_meshes_obs = Observable(Plot[])

        on(t) do _
            time_macro_step!(sim)

            for p in cell_meshes_obs[]
                delete!(ax1.scene, p)
            end
            cell_meshes_obs[] = Plot[]

            update_cell_colors!()

            cell_voxels = get_stat(sim.model, :pixels_by_cell)
            known_ids = Set(c.id for c in sim.cells)

            for (cid, voxels) in cell_voxels
                if !(cid in known_ids)
                    continue
                end
                cell = sim.cell_map[cid]

                color = get_cell_color(cell, sim) 
            
                p = draw_smooth_cell_mesh!(ax1, voxels, size(sim.model.grid), color)
                push!(cell_meshes_obs[], p)
            end       
            
            # Show selected vax only
            if show_fields && selected_vax[] != :none
                for (key, volplot) in volumes
                    volplot.visible[] = (key == selected_vax[])
                    if key == selected_vax[]
                        field = sim.vaxes[sim.cfg[:vaxes][key][:id]]
                        volplot[4][] = convert(Array{Float32, 3}, field.x)
                    end
                end
            end
            
            

            step_label.text[] = "Step: $(Int(round(t[])))"
            count_label.text[] = "Cell count: $(length(sim.cells))"
            
            pixels_by_cell = get_stat(sim.model, :pixels_by_cell)
            cell_volumes = [length(get(pixels_by_cell, c.id, [])) for c in sim.cells]
            
            vol_label.text[] = "Avg vol: $(round(mean(cell_volumes), digits=1))"
            max_vol_label.text[] = "Max vol: $(maximum(cell_volumes))"            
            
            # Count cells by cum_state
            state_counts = Dict{Symbol, Int}()
            for c in sim.cells
                state_sym = convert_idx_to_state(c.state[:cum_state], sim)
                state_counts[state_sym] = get(state_counts, state_sym, 0) + 1
            end

            state_text = join([string(k) * ": " * string(v) for (k, v) in sort(state_counts)], "   ")
            state_label.text[] = state_text
            
            


        end

    else
        fig = Figure()
        ax = fig[1, 1] = Axis(fig)
        vox = Observable(nothing)

        on(t) do _
            time_macro_step!(sim)
            vox[] = draw_cpm(sim,
                    num_cell_kinds = num_cell_kinds,
                    show_fields = show_fields,
                    show_activity = show_activity,
                    cell_kind_color = cell_kind_color,
                    field_color = field_color,
                    vox = vox[])
        end
    end

    return fig, t
end



"""
    anim_sim(sim::Sim, t::Observable; num_of_steps::Int = 1000, fps::Float64 = 1.0/30.0, should_we_continue = always_true)

Animates the simulation by incrementing time and redrawing plots interactively.

Arguments:
- `sim`: The simulation object
- `t`: A `Makie.Observable` that triggers frame updates (from `prep_vis_sim`)
- `num_of_steps`: Number of simulation steps to run
- `fps`: Time delay between frames (controls simulation speed)
- `should_we_continue`: Optional early stopping function `(sim) -> Bool`

Behavior:
- Increments `t[]` on each frame, triggering all associated visualization updates
- Sleeps for `fps` seconds between frames to create a real-time feel
- Used for live, exploratory simulations without recording

Returns nothing. Designed to be used inside `simulate`.
"""
function anim_sim(sim, t; num_of_steps = 1000, fps = 1.0/30., should_we_continue = always_true)
    for i in 1:num_of_steps
        t[] += 1
        if !should_we_continue(sim)
            break;
        end
        sleep(fps)
    end
end


"""
    record_sim(sim::Sim, t::Observable, img::SceneOrFigure;
               filename::String = "simulation.mp4",
               num_of_steps::Int = 100,
               fps::Float64 = 1.0/30.0,
               should_we_continue = always_true)

Records the simulation animation to an `.mp4` video file using GLMakie's `record` function.

Arguments:
- `sim`: The simulation object
- `t`: A `Makie.Observable` driving the animation
- `img`: The Makie figure returned by `prep_vis_sim`
- `filename`: Output filename for the video
- `num_of_steps`: Number of macrosteps (frames) to record
- `fps`: Frame rate (seconds per frame)
- `should_we_continue`: Optional early stopping condition

Behavior:
- Advances the simulation step-by-step while rendering each frame
- Updates `t[]` to trigger redraws
- Outputs an `.mp4` file with the given name

Intended for producing reproducible visualizations for presentations, reports, or supplementary videos.
"""
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

"""
    init_sim(cfg::Dict) -> Dict

Initializes and prepares a GBM spheroid simulation based on the provided configuration dictionary `cfg`.

This function:
- Creates the `Sim` object including chemical fields, CPM model, and cell metadata
- Builds initial cell population as defined in `cfg[:cells]`
- Sets up receptor arrays and internal mappings (`:vaxmap`, etc.)
- Initializes optional simulation elements like barriers
- Optionally prepares the interactive 2D or 3D GLMakie visualization (if `cfg[:runtime][:show_sim] == true`)

Returns a dictionary with:
- `:sim` — the simulation object
- `:img` — the Makie figure (if visualization is enabled)
- `:t` — a `Makie.Observable` time variable used for stepping the simulation interactively

This is the main entry point for launching a simulation via visual or headless execution.
"""

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

"""
    simulate(sim_desc::Dict; num_of_steps::Int = 1, record::Bool = false, filename::String = "simulation.mp4",
             fps::Float64 = 1.0/30.0, should_we_continue = always_true)

Runs the simulation described by `sim_desc`, which must be the result of `init_sim`.

Parameters:
- `num_of_steps`: Number of macrosteps to simulate.
- `record`: If true, records the simulation to a video file.
- `filename`: Output filename for the recording (only if `record = true`).
- `fps`: Frame rate for animation or recording.
- `should_we_continue`: Optional function `(sim) -> Bool` to allow early stopping.

Behavior:
- If visualization is enabled (`cfg[:runtime][:show_sim]`), the simulation is animated interactively or recorded.
- Otherwise, it runs headlessly without visualization.
"""


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

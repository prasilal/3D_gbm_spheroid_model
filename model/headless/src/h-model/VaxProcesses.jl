module VaxProcesses

const RREACT = 1
const RPROD = 2
const RK = 3
const RW = 4
const RR_ABSORB = 5

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

function perform_diffusion!(field, D, barrier_idxs)
    if D == 0 return end

    dst = zeros(Float64, size(field))
    M = fill(D, size(field))
    if barrier_idxs != nothing
        M[barrier_idxs] .= 0
    end
    if (field |> size |> length == 2)
        B = OffsetArray(masked_laplacian(2), -1:1,-1:1)
    else
        B = OffsetArray(masked_laplacian(3), -1:1,-1:1, -1:1)
    end
    masked_manifold_smoothen!(dst, field.x, M, B)
    field.x .+= dst
end

function perform_decay!(field, d)
    field.x .*= 1.0 - d
end

function correct_vals!(field)
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

function get_bound_vaxes(nof_vaxes, cells, cell_border_list, cells_border_pixels) # what receptors are on neighboring cells
    bound_vaxes = zeros(Float64, nof_vaxes)
    for id in (@>> cell_border_list keys collect filter(isneqv(BGID)))
        cidx = findfirst(c -> c.id == id, cells)
        if cidx != nothing
            bound_vaxes .+= cells[cidx].receptors .* (cell_border_list[id] / length(cells_border_pixels[id]))
        end
    end
    bound_vaxes
end

end # module

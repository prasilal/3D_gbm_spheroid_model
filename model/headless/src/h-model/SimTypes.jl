module SimTypes


"Background cell type index."
const BG_TYPE = BG_KIND
"Cell type representing dead cells."
const CELL_TYPE_DEAD = BG_TYPE + 1
"Starting index of active cell types."
const CELL_TYPES_START = CELL_TYPE_DEAD + 1
"Offset to add to cumulative state to get CPM kind index."
const CELL_TYPES_OFFSET = CELL_TYPE_DEAD

"Represents an individual cell in the simulation."
mutable struct Cell
    id :: Int
    receptors :: Vector{Float64}
    state :: Dict{Symbol, Any}
end

"Represents the overall simulation state."
mutable struct Sim{N}
    cells :: Vector{Cell}
    vaxes :: Vector{Grid{Float64, N}}
    model :: GridModel
    cfg :: Dict
    irep :: Dict{Symbol, Any}
end

end # module

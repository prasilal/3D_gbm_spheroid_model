module HModel

using Colors
const RGBAf0 = RGBA{Float32}
using GLMakie
using MLStyle
using Random
using Distributions
using Images
using GeometryBasics

include("../rf.jl")
include("../rf_cpm.jl")
include("../rf_cpm/Vizs.jl")
include("../rf_plot.jl")


include("SimTypes.jl")
include("SimBuilder.jl")
include("CPMParams.jl")
include("VaxProcesses.jl")
include("CellLogic.jl")
include("Rules.jl")
include("Visualization.jl")
include("SimControl.jl")
include("Util.jl")
include("InternalRep.jl")

using .SimTypes
using .SimBuilder
using .CPMParams
using .VaxProcesses
using .CellLogic
using .Rules
using .Visualization
using .SimControl
using .Util
using .InternalRep

end # module

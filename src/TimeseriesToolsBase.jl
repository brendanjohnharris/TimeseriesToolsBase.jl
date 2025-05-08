module TimeseriesToolsBase

import Unitful.unit

using Reexport
using DimensionalData
using IntervalSets
@reexport using DimensionalData
@reexport using IntervalSets

function __init__()
    ENV["UNITFUL_FANCY_EXPONENTS"] = true
end

include("Types.jl")
include("Utils.jl")
include("Spectra.jl")
include("Spectrograms.jl")
include("Unitful.jl")
include("Dates.jl")

end

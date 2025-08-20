module TimeseriesBase

using Reexport
using DimensionalData
using IntervalSets
@reexport using DimensionalData
@reexport using IntervalSets

include("ToolsArrays.jl")
@reexport using TimeseriesBase.ToolsArrays

include("TimeSeries.jl")
@reexport using TimeseriesBase.TimeSeries

include("Spectra.jl")
@reexport using TimeseriesBase.Spectra

include("UnitfulTools.jl")
@reexport using TimeseriesBase.UnitfulTools

include("Utils.jl")
@reexport using TimeseriesBase.Utils

include("DatesTools.jl")
@reexport using TimeseriesBase.DatesTools

include("Operators.jl")
@reexport using TimeseriesBase.Operators

end

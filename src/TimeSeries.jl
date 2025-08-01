module TimeSeries

export AbstractTimeSeries, AbstractTS,
       UnivariateTimeSeries, UnivariateTS,
       MultivariateTimeSeries, MultivariateTS,
       RegularTimeSeries, RegularTS,
       UnivariateRegular, MultivariateRegular,
       IrregularTimeSeries, IrregularTS,
       TimeIndex, RegularIndex, RegularTimeIndex,
       IrregularIndex, IrregularTimeIndex,
       TimeSeries, Timeseries,
       MultidimensionalIndex, MultidimensionalTimeSeries, MultidimensionalTS,
       SpikeTrain, MultivariateSpikeTrain, UnivariateSpikeTrain,
       spiketrain, spiketimes

using ..ToolsArrays
using DimensionalData
import DimensionalData: Dimension, TimeDim

Timeseries(x::DimArray) = ToolsArray(x)

"""
    TimeIndex

A type alias for a tuple containing a time dimension and any number of other dimensions.
"""
const TimeIndex = Tuple{A, Vararg{Dimension}} where {A <: TimeDim}

"""
    AbstractTimeSeries{T, N, B}

A type alias for an [AbstractDimArray](https://rafaqz.github.io/DimensionalData.jl/stable/api/#DimensionalData.AbstractDimArray) with a time index.
"""
const AbstractTimeSeries = AbstractTS = AbstractToolsArray{T, N, <:TimeIndex,
                                                           B} where {T, N, B}

"""
    UnivariateTimeSeries{T}

A type alias for a time series with one variable (a vector with only a `Ti` dimension).
"""
const UnivariateTimeSeries = UnivariateTS = AbstractTimeSeries{T, 1} where {T}

"""
    MultivariateTimeSeries{T}

A type alias for a multivariate time series (A matrix, with a first `Ti` dimension and an arbitrary second dimension).
"""
const MultivariateTimeSeries = MultivariateTS = AbstractTimeSeries{T, 2} where {T}

"""
    Var

A DimensionalData.jl dimension representing the variables of a multivariate time series.
"""
Var

"""
    RegularIndex

A type alias for a regularly sampled dimension, wrapping an `AbstractRange`.
"""
const RegularIndex = Dimensions.LookupArrays.Sampled{T, R} where {T, R <: AbstractRange}

"""
    RegularTimeIndex

A type alias for a tuple of dimensions containing a [`TimeIndex`](@ref) and any number of other dimensions.
"""
const RegularTimeIndex = Tuple{A,
                               Vararg{Dimension}} where {A <: TimeDim{<:RegularIndex}}

"""
    RegularTimeSeries{T, N, B}

A type alias for a regularly sampled time series.
"""
const RegularTimeSeries = RegularTS = AbstractToolsArray{T, N, <:RegularTimeIndex,
                                                         B} where {T, N, B}

const MultidimensionalIndex = Tuple{A,
                                    Vararg{Dimension{B}}} where {
                                                                 A <:
                                                                 TimeDim{<:RegularIndex},
                                                                 B <:
                                                                 RegularIndex
                                                                 }

"""
A multidimensional time series has a regular sampling over a dimension other than time; a one-dimensional time series can be thought of as a field over an even grid in 1 dimension that fluctuates over time.
"""
const MultidimensionalTimeSeries = AbstractToolsArray{T, N, <:MultidimensionalIndex,
                                                      B} where {T, N, B}
const MultidimensionalTS = MultidimensionalTimeSeries

"""
    IrregularIndex

A type alias for an irregularly sampled dimension, wrapping an `AbstractVector`.
"""
const IrregularIndex = Dimensions.LookupArrays.Sampled{T,
                                                       R} where {T,
                                                                 R <:
                                                                 AbstractVector
                                                                 }

"""
    IrregularTimeIndex

A type alias for a tuple of dimensions containing a [`TimeIndex`](@ref) and any number of other dimensions.
"""
const IrregularTimeIndex = Tuple{A,
                                 Vararg{Dimension}} where {A <:
                                                           TimeDim{<:IrregularIndex}}

"""
    IrregularTimeSeries

A type alias for a potentially irregularly sampled time series.
"""
const IrregularTimeSeries = IrregularTS = AbstractToolsArray{T, N, <:IrregularTimeIndex,
                                                             B} where {T, N, B}

"""
    BinaryTimeSeries

A type alias for a time series of bits.
"""
const BinaryTimeSeries = SpikeTrain = BinaryTS = AbstractToolsArray{T, N, <:TimeIndex,
                                                                    B} where {T <: Bool, N,
                                                                              B}

"""
    SpikeTrain

A type alias for a spike-train time series, which contains spike times in the time dimension and `true` for all values corresponding to a spike. The spike times can be retrieved with `times(x)`.
"""
SpikeTrain

const UnivariateSpikeTrain = typeintersect(UnivariateTimeSeries, SpikeTrain)
const MultivariateSpikeTrain = typeintersect(MultivariateTimeSeries, SpikeTrain)

function spiketrain(x; kwargs...)
    TimeSeries(sort(x), trues(length(x)); kwargs...)
end

function spiketimes(x::UnivariateSpikeTrain)
    times(x[x])
end
function spiketimes(x::SpikeTrain)
    map(spiketimes, eachslice(x, dims = tuple(2:ndims(x)...)))
end
spiketimes(x::AbstractArray) = x

"""
    Timeseries(x, t)

Constructs a univariate time series with time `t` and data `x`.

## Examples
```@example 1
julia> using TimeseriesToolsBase, Unitful;
julia> t = 1:100
julia> x = rand(100)
julia> ts = Timeseries(x, t)
julia> ts isa typeintersect(UnivariateTimeSeries, RegularTimeSeries)
```
"""
Timeseries(x, t; kwargs...) = ToolsArray(parent(x), (𝑡(t),); kwargs...)
Timeseries(x, t::TimeDim; kwargs...) = ToolsArray(parent(x), (t,); kwargs...)

"""
    TimeSeries(x, t, v)

Constructs a multivariate time series with time t, variable v, and data x.

## Examples
```@example 1
julia> t = 1:100;
julia> v = [:a, :b, :c];
julia> x = rand(100, 3);
julia> mts = Timeseries(x, t, v)
julia> mts isa typeintersect(MultivariateTimeSeries, RegularTimeSeries)
```
"""
Timeseries(x, t, v; kwargs...) = ToolsArray(x, (𝑡(t), Var(v)); kwargs...)

function Timeseries(x, t::TimeDim, dims::Vararg{<:Dimension}; kwargs...)
    ToolsArray(parent(x), (t, dims...); kwargs...)
end
function Timeseries(x, t, dims::Vararg{<:Dimension}; kwargs...)
    ToolsArray(parent(x), (𝑡(t), dims...); kwargs...)
end

convertconst(a, _) = a

const UnivariateRegular = typeintersect(UnivariateTimeSeries, RegularTimeSeries)
const MultivariateRegular = typeintersect(MultivariateTimeSeries, RegularTimeSeries)

end

using TestItems
using TestItemRunner

@run_package_tests

@testitem "Dates" begin
    using Dates, Unitful
    x = 1:100
    t = DateTime(1901):Year(1):DateTime(2000)
    y = @test_nowarn TimeSeries(t, x)
    @test y isa RegularTimeSeries
    @test samplingperiod(y) == Year(1)
    @test times(y) == t
    @test duration(y) == last(t) - first(t)
    @test unit(y) == NoUnits
end

@testitem "Spectra" begin
    using Unitful
    # Define a test time series
    fs = 0.1:0.1:100
    S = 1.0 ./ fs
    Pxx = Spectrum(fs, S)
    @test Pxx isa RegularSpectrum
    # ........and other funcs
end

@testitem "Unitful" begin
    using Unitful
    ts = (1:1000)u"s"
    x = @test_nowarn TimeSeries(ts, randn(1000))
    @test TimeSeries(ustripall(ts), collect(x), u"s") == x
    @test x isa AbstractTimeSeries
    @test x isa UnitfulTimeSeries
    @test x isa RegularTimeSeries
    @test x isa UnivariateTimeSeries

    @test step(x) == step(ts)
    @test samplingrate(x) == 1 / step(ts)
    @test times(x) == ts
    @test duration(x) == -first(-(extrema(ts)...))
    @test x[𝑡(1u"s" .. 10u"s")] == x[1:10]
    @test x[𝑡 = 1:10] == x[1:10]

    @test_nowarn rectify(x; dims = 𝑡)
end

include("Types.jl")
include("Utils.jl")
include("IO.jl")
include("Unitful.jl")
include("SpikeTrains.jl")
include("Operators.jl")
include("MakieExt.jl")
include("TimeseriesSurrogatesExt.jl")
include("Extensions.jl")

module Thermochron

    using Reexport
    @reexport using StatGeochemBase
    @reexport using NaNStatistics
    @reexport using Distributions

    using LinearAlgebra
    using LoopVectorization
    using ProgressMeter: Progress, update!, finish!
    using LsqFit: curve_fit
    using BangBang: setproperty!!

    const FloatRange = typeof(1.0:1.0:10.0)
    floatrange(x) = range(Float64(first(x)), Float64(last(x)), length=length(x))

    # Physical constants
    const SEC_MYR = 1E6*365.25*24*3600
    const LOG_SEC_MYR = log(SEC_MYR)

    # Decay constants
    const λ235U = log(2)/(7.0381*10^8)*10^6 # [1/Myr] Jaffey et al. 1971
    const λ238U = log(2)/(4.4683*10^9)*10^6 # [1/Myr] Jaffey et al. 1971
    const λ232Th = log(2)/(1.405*10^10)*10^6 # [1/Myr]
    const λ147Sm = log(2)/(1.070*10^11)*10^6 # [1/Myr] Kossert et al. 2009

    include("types.jl")
    include("utilities.jl")
    include("chronometers.jl")
    const ChronometerUnion{T} = Union{ZirconHe{T}, ApatiteHe{T}, ApatiteFT{T}, ApatiteTrackLength{T}}
    include("helium.jl")
    include("fissiontrack.jl")
    include("inversion.jl")
    include("show.jl")
end

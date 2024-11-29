module Thermochron

    using Reexport
    @reexport using StatGeochemBase
    @reexport using NaNStatistics
    @reexport using Distributions

    using LinearAlgebra
    using LoopVectorization
    using ProgressMeter: Progress, update!, finish!
    using LsqFit: curve_fit

    # Physical constants
    const SEC_MYR = 1E6*365.25*24*3600
    const LOG_SEC_MYR = log(SEC_MYR)

    # Jaffey decay constants
    const λ235U = log(2)/(7.0381*10^8)*10^6 # [1/Myr]
    const λ238U = log(2)/(4.4683*10^9)*10^6 # [1/Myr]
    const λ232Th = log(2)/(1.405*10^10)*10^6 # [1/Myr]

    include("types.jl")
    include("utilities.jl")
    include("chronometers.jl")
    include("helium.jl")
    include("fissiontrack.jl")
    include("inversion.jl")

end

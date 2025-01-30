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
    const λ235U = log(2)/(7.0381e8)*10^6    # [1/Myr] Jaffey et al. 1971
    const λ238U = log(2)/(4.4683e9)*10^6    # [1/Myr] Jaffey et al. 1971
    const λ232Th = log(2)/(1.405e10)*10^6   # [1/Myr] 
    const λ147Sm = log(2)/(1.070e11)*10^6   # [1/Myr] Kossert et al. 2009
    const λ40K = 5.5545e-10*10^6            # [1/Myr] Renne et al. 2010
    const λ40Kβ = 4.884e-10*10^6            # [1/Myr] Renne et al. 2010
    const λ40Kϵ = 0.580e-10*10^6            # [1/Myr] Renne et al. 2010
    const BR40K = λ40Kϵ/λ40K                # [unitless] branching ratio
    const κ40K = 1.1672e-4                  # [unitless] Garner et al., 1975 40K/K fraction
    export κ40K
    
    include("types.jl")
    export ZRDAAM, RDAAM                                                    # Damage-and-annealing based helium diffusivity model types
    export Ketcham1999FC, Ketcham2007FC, Yamada2005PC                       # Fission track annealing model types
    export Constraint, Unconformity, Boundary, DetailInterval               # Types used as inputs to MCMC functions

    include("chronometers.jl")
    const ChronometerUnion{T} = Union{ZirconFT{T}, ApatiteFT{T}, ApatiteTrackLength{T}, ZirconHe{T}, ApatiteHe{T}, GenericHe{T}, GenericAr{T}}
    export Chronometer, AbsoluteChronometer                             # Abstract types
    export ZirconFT, ApatiteFT, ApatiteTrackLength                      # Concrete fission track types
    export GenericHe, ZirconHe, ApatiteHe                               # Concrete U-Th/He types
    export GenericAr                                                    # Concrete K/Ar types
    export chronometers, empiricaluncertainty!, eU,                     # Functions
        get_age, get_age_sigma, set_age!, set_age_sigma!

    include("utilities.jl")
    include("helium.jl")
    include("argon.jl")
    export modelage                                                     # Functions

    include("fissiontrack.jl")
    export modellength                                                  # Functions

    include("inversion.jl")
    export MCMC, MCMC_varkinetics                                       # Functions

    include("show.jl")

    # Methodless functions for plotting extensions
    function ageeuplot end
    function ageeuplot! end
    export ageeuplot, ageeuplot!

end

module Thermochron

    using Reexport
    @reexport using Distributions
    @reexport using NaNStatistics
    @reexport using StatGeochemBase
    import StatGeochemBase.image_from_paths
    import StatGeochemBase.image_from_paths!

    using LinearAlgebra
    using LoopVectorization
    using ProgressMeter: Progress, update!, finish!
    using LsqFit: curve_fit
    using LogExpFunctions: logaddexp, logsubexp, logsumexp

    const FloatRange = typeof(1.0:1.0:10.0)
    floatrange(start, stop; length) = range(Float64(start), Float64(stop); length)
    floatrange(x::AbstractArray) = floatrange(first(x), last(x), length=length(x))

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
    export Constraint, Unconformity, Boundary, DetailInterval           # Types used as inputs to MCMC functions
    export Diffusivity
    
    include("chronometers.jl")
    const ChronometerUnion{T} = Union{ZirconFT{T}, MonaziteFT{T}, ApatiteFT{T}, ZirconTrackLength{T}, MonaziteTrackLength{T}, ApatiteTrackLength{T}, ZirconHe{T}, ApatiteHe{T}, SphericalHe{T}, PlanarHe{T}, SphericalAr{T}, PlanarAr{T}, MultipleDomain{T, SphericalAr{T}}, MultipleDomain{T, PlanarAr{T}}}
    export Chronometer, AbsoluteChronometer                             # Abstract types
    export ZirconTrackLength, MonaziteTrackLength, ApatiteTrackLength   # Concrete fission track length types
    export ZirconFT, MonaziteFT, ApatiteFT                              # Concrete fission track types
    export SphericalHe, PlanarHe, ZirconHe, ApatiteHe                   # Concrete U-Th/He types
    export SphericalAr, PlanarAr                                        # Concrete K/Ar types
    export chronometers, empiricaluncertainty!, eU,                     # Functions
        get_age, get_age_sigma, set_age!, set_age_sigma!

    include("argon.jl")
    include("helium.jl")
    export ZirconHeliumModel,ZRDAAM, ApatiteHeliumModel, RDAAM          # Damage-and-annealing based helium diffusivity model types

    include("mdd.jl")
    export MultipleDomain, MDDiffusivity

    include("fissiontrack.jl")
    export Ketcham1999FC, Ketcham2007FC                                 # Apatite fission track annealing model types
    export Yamada2007PC, Guenthner2013FC                                # Zircon fission track annealing models
    export Jones2021FA                                                  # Other mineral annealing models
    export modelage, modellength                                        # Functions
    const ModelUnion{T} = Union{Diffusivity{T}, MDDiffusivity{T}, RDAAM{T}, ZRDAAM{T}, Yamada2007PC{T}, Guenthner2013FC{T}, Ketcham1999FC{T}, Ketcham2007FC{T}, Jones2021FA{T}}

    include("utilities.jl")
    include("inversion.jl")
    export MCMC, MCMC_varkinetics                                       # Functions

    include("show.jl")

    # Methodless functions for plotting extensions
    function ageeuplot end
    function ageeuplot! end
    function errorbox end
    function errorbox! end
    export ageeuplot, ageeuplot!, errorbox, errorbox!

end

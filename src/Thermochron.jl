module Thermochron

    using Reexport
    @reexport using Distributions
    @reexport using NaNStatistics
    @reexport using StatGeochemBase
    import StatGeochemBase.image_from_paths
    import StatGeochemBase.image_from_paths!

    using Plots
    using Random
    using LinearAlgebra
    using ProgressMeter: Progress, update!, finish!
    using LsqFit: curve_fit
    using LogExpFunctions: logaddexp, logsubexp, logsumexp
    using LoopVectorization
    macro fast(ex)
        @static if Sys.isbsd()
            esc(:(@turbo check_empty=true $ex))
        else
            esc(:(@inbounds @simd ivdep $ex))
        end
    end

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
    
    # Type hierarchy, as well as some concrete types used within inversion functions
    include("types.jl")
    export Constraint, Boundary, DetailInterval, RegionalParameters     # Types used to pass information to MCMC functions
    
    # Chronometers objects used to represent different types of sample/measurement
    include("chronometers.jl")
    export Chronometer, AbsoluteChronometer                             # Abstract types
    export ZirconTrackLength, MonaziteTrackLength, ApatiteTrackLength, ApatiteTrackLengthOriented   # Concrete fission track length types
    export ZirconFT, MonaziteFT, ApatiteFT                              # Concrete fission track types
    export VitriniteRo                                                  # Concrete vitrinite reflectance type
    export SphericalHe, PlanarHe, ZirconHe, ApatiteHe, PlanarZirconHe   # Concrete U-Th-[Sm]/He types
    export SphericalAr, PlanarAr                                        # Concrete K/Ar types
    export get_age, get_age_sigma, empiricaluncertainty!, eU, radius    # Utility functions
    const PlanarNobleGas{T} = Union{PlanarZirconHe{T}, PlanarHe{T}, PlanarAr{T}}
    const SphericalNobleGas{T} = Union{ZirconHe{T}, ApatiteHe{T}, SphericalHe{T}, SphericalAr{T}}

    # Noble gas diffusion thermochronology functions and diffusivity/annealing models
    include("noblegas.jl")
    include("stepheating.jl")
    include("diffusion.jl")
    export Diffusivity, SDiffusivity                                    # Types for generic and scaled diffusivities
    export MDiffusivity, MSDiffusivity                                  # Types for multiple diffusivities
    export ZirconHeliumModel, ZRDAAM, ApatiteHeliumModel, RDAAM         # Damage-and-annealing based helium diffusivity model types
    export SingleDomain, MultipleDomain                                 # Types for modelling step-heating data, with one or more diffusion domains 
    export diffusivity                                                  # Function to calculate diffusivity at given temperature (and damage) for any DiffusivityModel

    # Fission track thermochronology functions and annealing models
    include("fissiontrack.jl")
    export Ketcham1999FC, Ketcham2007FC                                 # Apatite fission track annealing model types
    export Yamada2007PC, Guenthner2013FC                                # Zircon fission track annealing models
    export Jones2021FA                                                  # Other mineral annealing models
    export modelage, modellength, anneal!                               # Functions

    # Vitrinite reflectance thermochronology functions and annealing models
    include("vitrinite.jl")
    export NielsenBasinRo                                               # Vitrinite reflectance annealing models

    # Sample import
    include("parsing.jl")
    export chronometers, checktimediscretization                        # Parse datasets into Chronometer objects

    # Inversion
    include("utilities.jl")
    include("inversion.jl")
    include("show.jl")
    export MCMC, MCMC_varkinetics, model, model!                        # Inversion functions

    # Methodless functions for plotting extensions
    function ageeuplot end
    function ageeuplot! end
    function agesizeplot end
    function agesizeplot! end
    function errorbox end
    function errorbox! end
    export ageeuplot, ageeuplot!, agesizeplot, agesizeplot!, errorbox, errorbox!

end

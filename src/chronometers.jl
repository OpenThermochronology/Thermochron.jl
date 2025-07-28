# Abstract type to include any number of mineral chronometers (zircon, apatite, etc.)
abstract type Chronometer{T} end

# Implement methods to allow broadcasting
Base.length(x::Chronometer) = 1
Base.iterate(x::Chronometer) = (x, nothing)
Base.iterate(x::Chronometer, state) = nothing

# Implement methods to allow copying and comparison
Base.copy(x::Chronometer) = deepcopy(x)
Base.:(==)(x::Chronometer, y::Chronometer) = false
function Base.:(==)(x::T, y::T) where {T<:Chronometer}
    for n in fieldnames(T)
        isequal(getfield(x, n), getfield(y, n)) || return false
    end
    return true
end

# Abstract subtype for chronometers that include an absolute age and age uncertainty
abstract type AbsoluteChronometer{T} <:Chronometer{T} end  

# Abstract subtypes for different categories of chronometers
abstract type FissionTrackLength{T} <: Chronometer{T} end
abstract type FissionTrackSample{T} <: AbsoluteChronometer{T} end
abstract type HeliumSample{T} <: AbsoluteChronometer{T} end
abstract type ArgonSample{T} <: AbsoluteChronometer{T} end


## --- Fission track length types

"""
```julia
ApatiteTrackLengthOriented(T::Type{<:AbstractFloat}=Float64; 
    length::Number = NaN,                   # [um] fission track length
    angle::Number = NaN,                    # [degrees] track angle from the c-axis
    lcmod::Number = lcmod(length, angle),   # [um] model length of an equivalent c-axis parallel rack
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    l0::Number = 16.38,                     # [um] Initial track length
    l0_sigma::Number = 0.1311,              # [um] Initial track length unertainty
    dpar::Number = NaN,                     # [um] diameter parallel to track
    F::Number = NaN,                        # [APFU] F concentration, in atoms per formula unit
    Cl::Number = NaN,                       # [APFU] Cl concentration, in atoms per formula unit
    OH::Number = NaN,                       # [APFU] OH concentration, in atoms per formula unit
    rmr0::Number = NaN,                     # [unitless] annealing parameter
    ledges = (0:1.0:20),                    # [um] length bin edges, for internal model length histogram
    agesteps::AbstracVector | tsteps::AbstracVector, # Temporal discretization
)
```
Construct an `ApatiteTrackLengthOriented` chronometer representing a single apatite fission 
track `length` um long, oriented at `angle` degrees to the c-axis, with a relative annealing  
resistance specified by `rmr0`, optionally at a constant temperature offset (relative 
to other samples) of `offset` [C].

If not provided directly, `rmr0` will be calculated, in order of preference:
1. from `F`, `Cl`, and `OH` together, via the `rmr0model` function
2. from `Cl` alone, via the `rmr0fromcl` function
3. from `dpar`, via the `rmr0fromdpar` functions
4. using a default fallback value of 0.83, if none of the above are provided.

Temporal discretization follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` must be sorted 
in increasing order.
"""
struct ApatiteTrackLengthOriented{T<:AbstractFloat} <: FissionTrackLength{T}
    length::T               # [um] track length
    angle::T                # [degrees] track angle from the c-axis
    lcmod::T                # [um] model length of an equivalent c-axis parallel rack
    offset::T               # [C] temperature offset relative to other samples
    l0::T                   # [um] Initial track length
    l0_sigma::T             # [um] Initial track length unertainty
    agesteps::FloatRange    # [Ma] age in Ma relative to the present
    tsteps::FloatRange      # [Ma] forward time since crystallization
    r::Vector{T}            # [unitless]
    pr::Vector{T}           # [unitless]
    ledges::FloatRange      # [um] Length distribution edges
    ldist::Vector{T}        # [um] Length log likelihood
    rmr0::T                 # [unitless] relative resistance to annealing (0=most, 1=least)
end
function ApatiteTrackLengthOriented(T::Type{<:AbstractFloat}=Float64; 
        length::Number = NaN, 
        angle::Number = NaN, 
        lcmod::Number = lcmod(length, angle),
        offset::Number = zero(T),
        l0::Number = NaN,
        l0_sigma::Number = NaN,
        dpar::Number = NaN, 
        F::Number = NaN, 
        Cl::Number = NaN, 
        OH::Number = NaN, 
        rmr0::Number = NaN,
        ledges = (0:1.0:20),
        agesteps = nothing, 
        tsteps = nothing, 
    )
    # Temporal discretization
    tsteps, agesteps = checkdiscretization(tsteps, agesteps)
    # Multikinetic fission track parameters
    if isnan(rmr0)
        s = F + Cl + OH
        rmr0 = if !isnan(s)
            rmr0model(F/s*2, Cl/s*2, OH/s*2)
        elseif !isnan(Cl)
            rmr0fromcl(Cl)
        elseif !isnan(dpar)
            rmr0fromdpar(dpar)
        else
            0.83
        end
    end
    if isnan(l0) 
        l0 = if !isnan(dpar)
            apatitel0modfromdpar(dpar)
        else
            16.38
        end
    end
    if isnan(l0_sigma)
        l0_sigma = 0.1311
    end
    r=zeros(T, size(agesteps))
    pr=zeros(T, size(agesteps))
    ldist=zeros(T, size(ledges).-1)
    ApatiteTrackLengthOriented(
        T(length),
        T(angle),
        T(lcmod),
        T(offset),
        T(l0),
        T(l0_sigma),
        floatrange(agesteps),
        floatrange(tsteps),
        r,
        pr,
        floatrange(ledges),
        ldist,
        T(rmr0),
    )
end


"""
```julia
ApatiteTrackLength(T::Type{<:AbstractFloat}=Float64; 
    length::Number = NaN,                   # [um] fission track length
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    l0::Number = 16.38,                     # [um] Initial track length
    l0_sigma::Number = 0.1311,              # [um] Initial track length unertainty
    dpar::Number = NaN,                     # [um] diameter parallel to track
    F::Number = NaN,                        # [APFU] F concentration, in atoms per formula unit
    Cl::Number = NaN,                       # [APFU] Cl concentration, in atoms per formula unit
    OH::Number = NaN,                       # [APFU] OH concentration, in atoms per formula unit
    rmr0::Number = NaN,                     # [unitless] annealing parameter
    ledges = (0:1.0:20),                    # [um] length bin edges, for internal model length histogram
    agesteps::AbstracVector | tsteps::AbstracVector, # Temporal discretization
)
```
Construct an `ApatiteTrackLengthOriented` chronometer representing a single apatite 
fission track `length` um long with a relative annealing resistance specified by `rmr0`, 
optionally at a constant temperature offset (relative to other samples) of `offset` [C].

If not provided directly, `rmr0` will be calculated, in order of preference:
1. from `F`, `Cl`, and `OH` together, via the `rmr0model` function
2. from `Cl` alone, via the `rmr0fromcl` function
3. from `dpar`, via the `rmr0fromdpar` functions
4. using a default fallback value of 0.83, if none of the above are provided.

Temporal discretization follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` must be sorted 
in increasing order.
"""
struct ApatiteTrackLength{T<:AbstractFloat} <: FissionTrackLength{T}
    length::T               # [um] track length
    offset::T               # [C] temperature offset relative to other samples
    l0::T                   # [um] Initial track length
    l0_sigma::T             # [um] Initial track length unertainty
    agesteps::FloatRange    # [Ma] age in Ma relative to the present
    tsteps::FloatRange      # [Ma] forward time since crystallization
    r::Vector{T}            # [unitless]
    pr::Vector{T}           # [unitless]
    ledges::FloatRange      # [um] Length distribution edges
    ldist::Vector{T}        # [um] Length log likelihood
    rmr0::T                 # [unitless] relative resistance to annealing (0=most, 1=least)
end
function ApatiteTrackLength(T::Type{<:AbstractFloat}=Float64; 
        length::Number = NaN, 
        offset::Number = zero(T),
        l0::Number = NaN,
        l0_sigma::Number = NaN,
        dpar::Number = NaN, 
        F::Number = NaN, 
        Cl::Number = NaN, 
        OH::Number = NaN, 
        rmr0::Number = NaN,
        ledges = (0:1.0:20),
        agesteps = nothing, 
        tsteps = nothing, 
    )
    # Temporal discretization
    tsteps, agesteps = checkdiscretization(tsteps, agesteps)
    # Multikinetic fission track parameters
    if isnan(rmr0)
        s = F + Cl + OH
        rmr0 = if !isnan(s)
            rmr0model(F/s*2, Cl/s*2, OH/s*2)
        elseif !isnan(Cl)
            rmr0fromcl(Cl)
        elseif !isnan(dpar)
            rmr0fromdpar(dpar)
        else
            0.83
        end
    end
    if isnan(l0) 
        l0 = if !isnan(dpar)
            apatitel0modfromdpar(dpar)
        else
            16.38
        end
    end
    if isnan(l0_sigma)
        l0_sigma = 0.1311
    end
    r=zeros(T, size(agesteps))
    pr=zeros(T, size(agesteps))
    ldist=zeros(T, size(ledges).-1)
    ApatiteTrackLength(
        T(length),
        T(offset),
        T(l0),
        T(l0_sigma),
        floatrange(agesteps),
        floatrange(tsteps),
        r,
        pr,
        floatrange(ledges),
        ldist,
        T(rmr0),
    )
end


"""
```julia
ZirconTrackLength(T::Type{<:AbstractFloat}=Float64; 
    length::Number = NaN,                   # [um] fission track length
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    l0::Number = 11.17,                     # [um] Initial track length
    l0_sigma::Number = 0.051,               # [um] Initial track length unertainty    
    ledges = (0:1.0:20),                    # [um] length bin edges, for internal model length histogram
    agesteps::AbstracVector | tsteps::AbstracVector, # Temporal discretization
)
```
Construct a `ZirconTrackLength` chronometer representing a single zircon fission track
`length` um long, optionally at a constant temperature offset (relative to other  
samples) of `offset` [C].

Temporal discretization follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` must be sorted 
in increasing order.
"""
struct ZirconTrackLength{T<:AbstractFloat} <: FissionTrackLength{T}
    length::T               # [um] track length
    offset::T               # [C] temperature offset relative to other samples
    l0::T                   # [um] initial track length
    l0_sigma::T             # [um] initial track length uncertainty
    agesteps::FloatRange    # [Ma] age in Ma relative to the present
    tsteps::FloatRange      # [Ma] forward time since crystallization
    r::Vector{T}            # [unitless]
    pr::Vector{T}           # [unitless]
    ledges::FloatRange      # [um] Length distribution edges
    ldist::Vector{T}        # [um] Length log likelihood
end
function ZirconTrackLength(T::Type{<:AbstractFloat}=Float64; 
        length::Number = NaN, 
        offset::Number = zero(T),
        l0::Number = 11.17,
        l0_sigma::Number = 0.051,
        ledges = (0:1.0:20),
        agesteps = nothing, 
        tsteps = nothing,
    )
    # Temporal discretization
    tsteps, agesteps = checkdiscretization(tsteps, agesteps)
    # Initial track length and uncertainty
    if isnan(l0) 
        l0 = 11.17
    end
    if isnan(l0_sigma)
        l0_sigma = 0.051
    end
    r=zeros(T, size(agesteps))
    pr=zeros(T, size(agesteps))
    ldist=zeros(T, size(ledges).-1)
    ZirconTrackLength(
        T(length),
        T(offset),
        T(l0),
        T(l0_sigma),
        floatrange(agesteps),
        floatrange(tsteps),
        r,
        pr,
        floatrange(ledges),
        ldist,
    )
end


"""
```julia
MonaziteTrackLength(T::Type{<:AbstractFloat} = Float64; 
    length::Number = NaN,                   # [um] fission track length
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    l0::Number = 10.60,                     # [um] Initial track length
    l0_sigma::Number = 0.19,                # [um] Initial track length unertainty    
    ledges = (0:1.0:20),                    # [um] length bin edges, for internal model length histogram
    agesteps::AbstracVector | tsteps::AbstracVector, # Temporal discretization
)
```
Construct a `MonaziteTrackLength` chronometer representing a single monazite fission track
`length` um long, optionally at a constant temperature offset (relative to other  
samples) of `offset` [C].

Temporal discretization follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` must be sorted 
in increasing order.
"""
struct MonaziteTrackLength{T<:AbstractFloat} <: FissionTrackLength{T}
    length::T               # [um] track length
    offset::T               # [C] temperature offset relative to other samples
    l0::T                   # [um] initial track length
    l0_sigma::T             # [um] initial track length uncertainty
    agesteps::FloatRange    # [Ma] age in Ma relative to the present
    tsteps::FloatRange      # [Ma] forward time since crystallization
    r::Vector{T}            # [unitless]
    pr::Vector{T}           # [unitless]
    ledges::FloatRange      # [um] Length distribution edges
    ldist::Vector{T}        # [um] Length log likelihood
end
function MonaziteTrackLength(T::Type{<:AbstractFloat}=Float64; 
        length::Number = NaN,
        offset::Number = zero(T),
        l0::Number = 10.60,
        l0_sigma::Number = 0.19,
        ledges = (0:1.0:20),
        agesteps = nothing, 
        tsteps = nothing, 
    )
    # Temporal discretization
    tsteps, agesteps = checkdiscretization(tsteps, agesteps)
    # Initial track length and uncertainty
    if isnan(l0) 
        l0 = 10.60
    end
    if isnan(l0_sigma)
        l0_sigma = 0.19
    end
    r=zeros(T, size(agesteps))
    pr=zeros(T, size(agesteps))
    ldist=zeros(T, size(ledges).-1)
    MonaziteTrackLength(
        T(length),
        T(offset),
        T(l0),
        T(l0_sigma),
        floatrange(agesteps),
        floatrange(tsteps),
        r,
        pr,
        floatrange(ledges),
        ldist,
    )
end

## --- Fission track age types

"""
```julia
ZirconFT(T::Type{<:AbstractFloat} = Float64; 
    age::Number = NaN,              # [Ma] fission track age
    age_sigma::Number = NaN,        # [Ma] fission track age uncertainty
    offset::Number = zero(T),       # [C] temperature offset relative to other samples
    agesteps::AbstracVector | tsteps::AbstracVector, 
)
```
Construct a `ZirconFT` chronometer representing a zircon with a fission track age 
of `age` ± `age_sigma` [Ma], optionally at a constant temperature offset (relative to 
other samples) of `offset` [C].

Temporal discretization follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` must be sorted 
in increasing order.
"""
struct ZirconFT{T<:AbstractFloat} <: FissionTrackSample{T}
    age::T                  # [Ma] fission track age
    age_sigma::T            # [Ma] fission track age uncertainty (one-sigma)
    offset::T               # [C] temperature offset relative to other samples
    agesteps::FloatRange    # [Ma] age in Ma relative to the present
    tsteps::FloatRange      # [Ma] forward time since crystallization
end
function ZirconFT(T::Type{<:AbstractFloat}=Float64; 
        age::Number = NaN,              # [Ma] fission track age
        age_sigma::Number = NaN,        # [Ma] fission track age uncertainty
        offset::Number = zero(T),       # [C] temperature offset relative to other samples
        agesteps = nothing, 
        tsteps = nothing, 
    )
    # Temporal discretization
    tsteps, agesteps = checkdiscretization(tsteps, agesteps)
    ZirconFT(
        T(age),
        T(age_sigma),
        T(offset),
        floatrange(agesteps),
        floatrange(tsteps),
    )
end


"""
```julia
MonaziteFT(T::Type{<:AbstractFloat} = Float64; 
    age::Number = NaN,              # [Ma] fission track age
    age_sigma::Number = NaN,        # [Ma] fission track age uncertainty
    offset::Number = zero(T),       # [C] temperature offset relative to other samples
    agesteps::AbstracVector | tsteps::AbstracVector, # Temporal discretization
)
```
Construct a `MonaziteFT` chronometer representing a monazite with a fission track age 
of `age` ± `age_sigma` [Ma], optionally at a constant temperature offset (relative to 
other samples) of `offset` [C].

Temporal discretization follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` must be sorted 
in increasing order.
"""
struct MonaziteFT{T<:AbstractFloat} <: FissionTrackSample{T}
    age::T                  # [Ma] fission track age
    age_sigma::T            # [Ma] fission track age uncertainty (one-sigma)
    offset::T               # [C] temperature offset relative to other samples
    agesteps::FloatRange    # [Ma] age in Ma relative to the present
    tsteps::FloatRange      # [Ma] forward time since crystallization
end
function MonaziteFT(T::Type{<:AbstractFloat}=Float64; 
        age::Number = NaN,              # [Ma] fission track age
        age_sigma::Number = NaN,        # [Ma] fission track age uncertainty
        offset::Number = zero(T),       # [C] temperature offset relative to other samples
        agesteps=nothing,
        tsteps=nothing, 
    )
    # Temporal discretization
    tsteps, agesteps = checkdiscretization(tsteps, agesteps)
    MonaziteFT(
        T(age),
        T(age_sigma),
        T(offset),
        floatrange(agesteps),
        floatrange(tsteps),
    )
end


"""
```julia
ApatiteFT(T::Type{<:AbstractFloat} = Float64; 
    age::Number = NaN,              # [Ma] fission track age
    age_sigma::Number = NaN,        # [Ma] fission track age uncertainty
    offset::Number = zero(T),       # [C] temperature offset relative to other samples
    dpar::Number = NaN,             # [um] diameter parallel to track
    F::Number = NaN,                # [APFU] F concentration, in atoms per formula unit
    Cl::Number = NaN,               # [APFU] Cl concentration, in atoms per formula unit
    OH::Number = NaN,               # [APFU] OH concentration, in atoms per formula unit
    rmr0::Number = NaN,             # [unitless] annealing parameter
    agesteps::AbstracVector | tsteps::AbstracVector, # Temporal discretization
)
```
Construct an `ApatiteFT` chronometer representing a apatite with a fission track age 
of `age` ± `age_sigma` [Ma] and a relative annealing resistance specified by `rmr0`,
and optionally a constant temperature offset (relative to other samples) of `offset` [C].

If not provided directly, `rmr0` will be calculated, in order of preference:
1. from `F`, `Cl`, and `OH` together, via the `rmr0model` function
2. from `Cl` alone, via the `rmr0fromcl` function
3. from `dpar`, via the `rmr0fromdpar` functions
4. using a default fallback value of 0.83, if none of the above are provided.

Temporal discretization follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` must be sorted 
in increasing order.
"""
struct ApatiteFT{T<:AbstractFloat} <: FissionTrackSample{T}
    age::T                  # [Ma] fission track age
    age_sigma::T            # [Ma] fission track age uncertainty (one-sigma)
    offset::T               # [C] temperature offset relative to other samples
    agesteps::FloatRange    # [Ma] age in Ma relative to the present
    tsteps::FloatRange      # [Ma] forward time since crystallization
    rmr0::T                 # [unitless] relative resistance to annealing (0=most, 1=least)
end
function ApatiteFT(T::Type{<:AbstractFloat}=Float64; 
        age::Number = NaN, 
        age_sigma::Number = NaN, 
        offset::Number = zero(T),
        dpar::Number = NaN,
        F::Number = NaN,
        Cl::Number = NaN,
        OH::Number = NaN,
        rmr0::Number = NaN,
        agesteps=nothing, 
        tsteps=nothing, 
    )
    # Temporal discretization
    tsteps, agesteps = checkdiscretization(tsteps, agesteps)
    # Multikinetic fission track parameters
    if isnan(rmr0)
        s = F + Cl + OH
        rmr0 = if !isnan(s)
            rmr0model(F/s*2, Cl/s*2, OH/s*2)
        elseif !isnan(Cl)
            rmr0fromcl(Cl)
        elseif !isnan(dpar)
            rmr0fromdpar(dpar)
        else
            0.83
        end
    end
    ApatiteFT(
        T(age),
        T(age_sigma),
        T(offset),
        floatrange(agesteps),
        floatrange(tsteps),
        T(rmr0),
    )
end

## --- Helium sample types

"""
```julia
ZirconHe(T=Float64;
    age::Number = T(NaN),                   # [Ma] raw helium age
    age_sigma::Number = T(NaN),             # [Ma] raw helium age uncertainty (one-sigma)
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    r::Number,                              # [um] spherical radius 
    dr::Number = one(T),                    # [um] radial step size
    U238::Number,                           # [ppm] zircon U-238 concentration
    Th232::Number,                          # [ppm] zircon Th-232 concentration
    Sm147::Number = zero(T),                # [ppm] zircon Sm-147 concentration
    U238_matrix::Number = zero(T),          # [ppm] matrix U-238 concentration
    Th232_matrix::Number = zero(T),         # [ppm] matrix Th-232 concentration
    Sm147_matrix::Number = zero(T),         # [ppm] matrix Sm-147 concentration
    volumeweighting::Symbol=:cylindrical,   # (:spherical, :cylindrical, or :planar) relative volume proportions of each radial model shell, for averaging purposes
    agesteps::AbstractVector | tsteps::AbstractVector, # Temporal discretization
)
```
Construct a `ZirconHe` chronometer representing a zircon with a raw 
helium age of `age` ± `age_sigma` [Ma], a spherical radius of `r` [μm], and  
uniform U, Th and Sm concentrations specified by `U238`, `Th232`, and `Sm147` [PPMw]. 
A present day U-235/U-238 ratio of 1/137.818 is assumed.

Spatial discretization follows a halfwidth step of `dr` [μm], while temporal 
discretization  follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` 
must be sorted in increasing order.
"""
struct ZirconHe{T<:AbstractFloat} <: HeliumSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to other samples
    agesteps::FloatRange        # [Ma] age in Ma relative to the present
    tsteps::FloatRange          # [Ma] forward time since crystallization
    rsteps::FloatRange          # [um] radius bin centers
    redges::FloatRange          # [um] radius bin edges
    relvolumes::Vector{T}       # [unitless] fraction of volume in each radial step
    nrsteps::Int                # [n] number of radial steps, including both implicit points at each side
    r238U::Vector{T}            # [atoms/g] radial U-238 concentrations
    r235U::Vector{T}            # [atoms/g] radial U-235 concentrations
    r232Th::Vector{T}           # [atoms/g] radial Th-232 concentrations
    r147Sm::Vector{T}           # [atoms/g] radial Sm-147 concentrations
    alphadeposition::Matrix{T}  # [atoms/g] alpha deposition matrix
    alphadamage::Matrix{T}      # [decays/g] initial damage matrix
    pr::Matrix{T}               # [unitless] reduced damage density matrix
    annealeddamage::Matrix{T}   # [decays/g] annealed damage matrix
    u::Matrix{T}
    β::Vector{T}
    Dz::Vector{T}
    DN17::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{T, Tridiagonal{T, Vector{T}}, Vector{Int64}}
    y::Vector{T}
end
function ZirconHe(T::Type{<:AbstractFloat}=Float64;
        age::Number = T(NaN),
        age_sigma::Number = T(NaN),
        offset::Number = zero(T),
        r::Number,
        dr::Number = one(T), 
        U238::Number,
        Th232::Number,
        Sm147::Number = zero(T),
        U238_matrix::Number = zero(T), 
        Th232_matrix::Number = zero(T), 
        Sm147_matrix::Number = zero(T), 
        volumeweighting::Symbol=:cylindrical,
        agesteps = nothing,
        tsteps = nothing,
    )

    # Temporal discretization
    tsteps, agesteps = checkdiscretization(tsteps, agesteps)
    tsteps, agesteps = floatrange(tsteps), floatrange(agesteps)

    # Zircon alpha stopping distances for each isotope in each decay chain, from
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alpharadii238U = (11.78, 14.09, 13.73, 14.13, 17.32, 16.69, 28.56, 16.48,)
    alpharadii235U = (12.58, 15.04, 19.36, 18.06, 23.07, 26.87, 22.47,)
    alpharadii232Th = (10.99, 16.67, 18.16, 17.32, 23.61, 29.19,)
    # Ketcham et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alpharadii147Sm = (4.76,)

    # Crystal size and spatial discretization
    redges = floatrange(0 : dr : r)                 # Edges of each radius element
    rsteps = cntr(redges)                           # Centers of each radius element
    nrsteps = length(rsteps)+2                      # Number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    # Relative volume each concentric radial shell
    relvolumes = volumefraction(volumeweighting, redges, r)

    # Additional discretization outside of grain, for alpha injection
    outsideredges = floatrange(r : dr : r+maximum(alpharadii238U))
    outsidersteps = cntr(outsideredges)
    outsiderelvolumes = volumefraction(volumeweighting, outsideredges, r)

    # Observed radial HPE profiles at present day
    r238U = fill(T(U238), size(rsteps))         # [PPMw]
    r235U = fill(T(U238/137.818), size(rsteps)) # [PPMw]
    r232Th = fill(T(Th232), size(rsteps))       # [PPMw]
    r147Sm = fill(T(Sm147), size(rsteps))       # [PPMw]

    # Convert to atoms per gram
    r238U .*= 6.022E23 / 1E6 / 238
    r235U .*= 6.022E23 / 1E6 / 235
    r232Th .*= 6.022E23 / 1E6 / 232
    r147Sm .*= 6.022E23 / 1E6 / 147

    # Outside (bulk/matrix) HPE concentrations, in atoms per gram
    o238U = U238_matrix * 6.022E23 / 1E6 / 238
    o235U = U238_matrix/137.818 * 6.022E23 / 1E6 / 235
    o232Th = Th232_matrix * 6.022E23 / 1E6 / 232
    o147Sm = Sm147_matrix * 6.022E23 / 1E6 / 147

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance

    # Allocate intersection density vector
    dint = zeros(T, length(redges) - 1)
    
    # Radial alpha deposition from U-238
    r238UHe = zeros(T, size(r238U))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, relvolumes, r238U)
        for alpharadius in alpharadii238U
            sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,rsteps[ri])
            @. r238UHe += relvolumes[ri] * dint * r238U[ri]
        end
    end
    # Correct for alpha injection
    if o238U > 0
        @inbounds for ri in eachindex(outsidersteps, outsiderelvolumes)
            for alpharadius in alpharadii238U
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,outsidersteps[ri])
                @. r238UHe += outsiderelvolumes[ri] * dint * o238U
            end
        end
    end

    # Radial alpha deposition from U-235
    r235UHe = zeros(T, size(r235U))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, relvolumes, r235U)
        for alpharadius in alpharadii235U
            sphereintersectiondensity!(dint, redges,relvolumes,alpharadius,rsteps[ri])
            @. r235UHe += relvolumes[ri] * dint * r235U[ri]
        end
    end
    # Correct for alpha injection
    if o235U > 0
        @inbounds for ri in eachindex(outsidersteps, outsiderelvolumes)
            for alpharadius in alpharadii235U
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,outsidersteps[ri])
                @. r235UHe += outsiderelvolumes[ri] * dint * o235U
            end
        end
    end

    # Radial alpha deposition from Th-232
    r232ThHe = zeros(T, size(r232Th))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, relvolumes, r232Th)
        for alpharadius in alpharadii232Th
            sphereintersectiondensity!(dint, redges,relvolumes,alpharadius,rsteps[ri])
            @. r232ThHe += relvolumes[ri] * dint * r232Th[ri]
        end
    end
    # Correct for alpha injection
    if o232Th > 0
        @inbounds for ri in eachindex(outsidersteps, outsiderelvolumes)
            for alpharadius in alpharadii232Th
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,outsidersteps[ri])
                @. r232ThHe += outsiderelvolumes[ri] * dint * o232Th
            end
        end
    end

    # Radial alpha deposition from Sm-147
    r147SmHe = zeros(T, size(r147Sm))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, relvolumes, r147Sm)
        for alpharadius in alpharadii147Sm
            sphereintersectiondensity!(dint, redges,relvolumes,alpharadius,rsteps[ri])
            @. r147SmHe += relvolumes[ri] * dint * r147Sm[ri]
        end
    end
    # Correct for alpha injection
    if o147Sm > 0
        @inbounds for ri in eachindex(outsidersteps, outsiderelvolumes)
            for alpharadius in alpharadii147Sm
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,outsidersteps[ri])
                @. r147SmHe += outsiderelvolumes[ri] * dint * o147Sm
            end
        end
    end   

    # Alpha decay recoil damage
    r238Udam = 8*r238U # No smoothing of alpha damage, 8 alphas per 238U
    r235Udam = 7*r235U # No smoothing of alpha damage, 7 alphas per 235U
    r232Thdam = 6*r232Th # No smoothing of alpha damage, 6 alphas per 232 Th
    r147Smdam = 1*r147Sm # No smoothing of alpha damage, 1 alpha per 147 Sm

    # Calculate corrected alpha deposition and recoil damage each time step for each radius
    dt_2 = step(tsteps)/2
    decay = zeros(T, length(tsteps))
    # Allocate deposition and damage arrays
    alphadeposition = zeros(T, length(tsteps), nrsteps-2)
    alphadamage = zeros(T, length(tsteps), nrsteps-2)
    pr = zeros(T, length(tsteps), length(tsteps))

    # U-238
    @. decay = exp(λ238U*(agesteps + dt_2)) - exp(λ238U*(agesteps - dt_2))
    mul!(alphadeposition, decay, r238UHe', one(T), one(T))
    mul!(alphadamage, decay, r238Udam', one(T), one(T))
    # U-235
    @. decay = exp(λ235U*(agesteps + dt_2)) - exp(λ235U*(agesteps - dt_2))
    mul!(alphadeposition, decay, r235UHe', one(T), one(T))
    mul!(alphadamage, decay, r235Udam', one(T), one(T))
    # Th-232
    @. decay = exp(λ232Th*(agesteps + dt_2)) - exp(λ232Th*(agesteps - dt_2))
    mul!(alphadeposition, decay, r232ThHe', one(T), one(T))
    mul!(alphadamage, decay, r232Thdam', one(T), one(T))
    # Sm-147
    @. decay = exp(λ147Sm*(agesteps + dt_2)) - exp(λ147Sm*(agesteps - dt_2))
    mul!(alphadeposition, decay, r147SmHe', one(T), one(T))
    mul!(alphadamage, decay, r147Smdam', one(T), one(T))

    # Allocate additional variables that will be needed for Crank-Nicolson
    annealeddamage = similar(alphadamage)
    β = zeros(T, nrsteps)

    # Allocate arrays for diffusivities
    Dz = zeros(T, length(tsteps))
    DN17 = zeros(T, length(tsteps))

    # Allocate output matrix for all timesteps
    u = zeros(T, nrsteps, length(tsteps)+1)

    # Allocate variables for tridiagonal matrix and RHS
    dl = ones(T, nrsteps-1)    # Sub-diagonal row
    d = ones(T, nrsteps)       # Diagonal
    du = ones(T, nrsteps-1)    # Supra-diagonal row
    du2 = ones(T, nrsteps-2)   # sup-sup-diagonal row for pivoting

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = Tridiagonal(dl, d, du, du2)
    F = lu(A, allowsingular=true)

    # Vector for RHS of Crank-Nicolson equation with regular grid cells
    y = zeros(T, nrsteps)

    return ZirconHe(
        T(age),
        T(age_sigma),
        T(offset),
        agesteps,
        tsteps,
        rsteps,
        redges,
        relvolumes,
        nrsteps,
        r238U,
        r235U,
        r232Th,
        r147Sm,
        alphadeposition,
        alphadamage,
        pr,
        annealeddamage,
        u,
        β,
        Dz,
        DN17,
        A,
        F,
        y,
    )
end


"""
```julia
ApatiteHe(T=Float64;
    age::Number = T(NaN),                   # [Ma] raw helium age
    age_sigma::Number = T(NaN),             # [Ma] raw helium age uncertainty (one-sigma)
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    r::Number,                              # [um] spherical radius 
    dr::Number = one(T),                    # [um] radial step size
    U238::Number,                           # [ppm] apatite U-238 concentration
    Th232::Number,                          # [ppm] apatite Th-232 concentration
    Sm147::Number = zero(T),                # [ppm] apatite Sm-147 concentration
    U238_matrix::Number = zero(T),          # [ppm] matrix U-238 concentration
    Th232_matrix::Number = zero(T),         # [ppm] matrix Th-232 concentration
    Sm147_matrix::Number = zero(T),         # [ppm] matrix Sm-147 concentration
    volumeweighting::Symbol=:cylindrical,   # (:spherical, :cylindrical, or :planar) relative volume proportions of each radial model shell, for averaging purposes    agesteps::AbstractVector | tsteps::AbstractVector, # Temporal discretization
)
```
Construct an `ApatiteHe` chronometer representing an apatite with a raw 
helium age of `age` ± `age_sigma` [Ma], a spherical radius of `r` [μm], and  
uniform U, Th and Sm concentrations specified by `U238`, `Th232`, and `Sm147` [PPMw]. 
A present day U-235/U-238 ratio of 1/137.818 is assumed.

Spatial discretization follows a halfwidth step of `dr` [μm], while temporal 
discretization  follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` 
must be sorted in increasing order.
"""
struct ApatiteHe{T<:AbstractFloat} <: HeliumSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to other samples
    agesteps::FloatRange        # [Ma] age in Ma relative to the present
    tsteps::FloatRange          # [Ma] forward time since crystallization
    rsteps::FloatRange          # [um] radius bin centers
    redges::FloatRange          # [um] radius bin edges
    relvolumes::Vector{T}       # [unitless] fraction of volume in each radial step
    nrsteps::Int                # [n] number of radial steps, including both implicit points at each side
    r238U::Vector{T}            # [atoms/g] radial U-238 concentrations
    r235U::Vector{T}            # [atoms/g] radial U-235 concentrations
    r232Th::Vector{T}           # [atoms/g] radial Th-232 concentrations
    r147Sm::Vector{T}           # [atoms/g] radial Sm-147 concentrations
    alphadeposition::Matrix{T}  # [atoms/g] alpha deposition matrix
    alphadamage::Matrix{T}      # [decays/g] initial damage matrix
    pr::Matrix{T}               # [unitless] reduced damage density matrix
    annealeddamage::Matrix{T}   # [decays/g] annealed damage matrix
    u::Matrix{T}
    β::Vector{T}
    DL::Vector{T}
    Dtrap::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{T, Tridiagonal{T, Vector{T}}, Vector{Int64}}
    y::Vector{T}
end
function ApatiteHe(T::Type{<:AbstractFloat}=Float64;
        age::Number = T(NaN),
        age_sigma::Number = T(NaN),
        offset::Number = zero(T),
        r::Number, 
        dr::Number = one(T), 
        U238::Number, 
        Th232::Number, 
        Sm147::Number = zero(T), 
        U238_matrix::Number = zero(T), 
        Th232_matrix::Number = zero(T), 
        Sm147_matrix::Number = zero(T), 
        volumeweighting::Symbol = :cylindrical,
        agesteps = nothing,
        tsteps = nothing,
    )

    # Temporal discretization
    tsteps, agesteps = checkdiscretization(tsteps, agesteps)
    tsteps, agesteps = floatrange(tsteps), floatrange(agesteps)

    # Apatite alpha stopping distances for each isotope in each decay chain, from
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alpharadii238U = (13.54, 16.26, 15.84, 16.31, 20.09, 22.89, 33.39, 19.10,)
    alpharadii235U = (14.48, 17.39, 22.5, 20.97, 26.89, 31.40, 26.18,)
    alpharadii232Th = (12.60, 19.32, 21.08, 20.09, 27.53, 34.14,)
    # Ketcham et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alpharadii147Sm = (5.93,)

    # Crystal size and spatial discretization
    redges = floatrange(0 : dr : r)                 # Edges of each radius element
    rsteps = cntr(redges)                           # Centers of each radius element
    nrsteps = length(rsteps)+2                      # Number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    # Relative volume each concentric radial shell
    relvolumes = volumefraction(volumeweighting, redges, r)

    # Additional discretization outside of grain, for alpha injection
    outsideredges = floatrange(r : dr : r+maximum(alpharadii238U))
    outsidersteps = cntr(outsideredges)
    outsiderelvolumes = volumefraction(volumeweighting, outsideredges, r)

    # Observed radial HPE profiles at present day
    r238U = fill(T(U238), size(rsteps))         # [PPMw]
    r235U = fill(T(U238/137.818), size(rsteps)) # [PPMw]
    r232Th = fill(T(Th232), size(rsteps))       # [PPMw]
    r147Sm = fill(T(Sm147), size(rsteps))       # [PPMw]

    # Convert to atoms per gram
    r238U .*= 6.022E23 / 1E6 / 238
    r235U .*= 6.022E23 / 1E6 / 235
    r232Th .*= 6.022E23 / 1E6 / 232
    r147Sm .*= 6.022E23 / 1E6 / 147

    # Outside (bulk/matrix) HPE concentrations, in atoms per gram
    o238U = U238_matrix * 6.022E23 / 1E6 / 238
    o235U = U238_matrix/137.818 * 6.022E23 / 1E6 / 235
    o232Th = Th232_matrix * 6.022E23 / 1E6 / 232
    o147Sm = Sm147_matrix * 6.022E23 / 1E6 / 147

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance

    # Allocate intersection density vector
    dint = zeros(T, length(redges) - 1)

    # Radial alpha deposition from U-238
    r238UHe = zeros(T, size(r238U))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, relvolumes, r238U)
        for alpharadius in alpharadii238U
            sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,rsteps[ri])
            @. r238UHe += relvolumes[ri] * dint * r238U[ri]
        end
    end
    # Correct for alpha injection
    if o238U > 0
        @inbounds for ri in eachindex(outsidersteps, outsiderelvolumes)
            for alpharadius in alpharadii238U
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,outsidersteps[ri])
                @. r238UHe += outsiderelvolumes[ri] * dint * o238U
            end
        end
    end

    # Radial alpha deposition from U-235
    r235UHe = zeros(T, size(r235U))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, relvolumes, r235U)
        for alpharadius in alpharadii235U
            sphereintersectiondensity!(dint, redges,relvolumes,alpharadius,rsteps[ri])
            @. r235UHe += relvolumes[ri] * dint * r235U[ri]
        end
    end
    # Correct for alpha injection
    if o235U > 0
        @inbounds for ri in eachindex(outsidersteps, outsiderelvolumes)
            for alpharadius in alpharadii235U
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,outsidersteps[ri])
                @. r235UHe += outsiderelvolumes[ri] * dint * o235U
            end
        end
    end

    # Radial alpha deposition from Th-232
    r232ThHe = zeros(T, size(r232Th))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, relvolumes, r232Th)
        for alpharadius in alpharadii232Th
            sphereintersectiondensity!(dint, redges,relvolumes,alpharadius,rsteps[ri])
            @. r232ThHe += relvolumes[ri] * dint * r232Th[ri]
        end
    end
    # Correct for alpha injection
    if o232Th > 0
        @inbounds for ri in eachindex(outsidersteps, outsiderelvolumes)
            for alpharadius in alpharadii232Th
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,outsidersteps[ri])
                @. r232ThHe += outsiderelvolumes[ri] * dint * o232Th
            end
        end
    end

    # Radial alpha deposition from Sm-147
    r147SmHe = zeros(T, size(r147Sm))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, relvolumes, r147Sm)
        for alpharadius in alpharadii147Sm
            sphereintersectiondensity!(dint, redges,relvolumes,alpharadius,rsteps[ri])
            @. r147SmHe += relvolumes[ri] * dint * r147Sm[ri]
        end
    end
    # Correct for alpha injection
    if o147Sm > 0
        @inbounds for ri in eachindex(outsidersteps, outsiderelvolumes)
            for alpharadius in alpharadii147Sm
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,outsidersteps[ri])
                @. r147SmHe += outsiderelvolumes[ri] * dint * o147Sm
            end
        end
    end

    # Alpha decay recoil damage
    r238Udam = 8*r238U # No smoothing of alpha damage, 8 alphas per 238U
    r235Udam = 7*r235U # No smoothing of alpha damage, 7 alphas per 235U
    r232Thdam = 6*r232Th # No smoothing of alpha damage, 6 alphas per 232 Th
    r147Smdam = 1*r147Sm # No smoothing of alpha damage, 1 alpha per 147 Sm

    # Calculate corrected alpha deposition and recoil damage each time step for each radius
    dt_2 = step(tsteps)/2
    decay = zeros(T, length(tsteps))
    # Allocate deposition and damage arrays
    alphadeposition = zeros(T, length(tsteps), nrsteps-2)
    alphadamage = zeros(T, length(tsteps), nrsteps-2)
    pr = zeros(T, length(tsteps), length(tsteps))

    # U-238
    @. decay = exp(λ238U*(agesteps + dt_2)) - exp(λ238U*(agesteps - dt_2))
    mul!(alphadeposition, decay, r238UHe', one(T), one(T))
    mul!(alphadamage, decay, r238Udam', one(T), one(T))
    # U-235
    @. decay = exp(λ235U*(agesteps + dt_2)) - exp(λ235U*(agesteps - dt_2))
    mul!(alphadeposition, decay, r235UHe', one(T), one(T))
    mul!(alphadamage, decay, r235Udam', one(T), one(T))
    # Th-232
    @. decay = exp(λ232Th*(agesteps + dt_2)) - exp(λ232Th*(agesteps - dt_2))
    mul!(alphadeposition, decay, r232ThHe', one(T), one(T))
    mul!(alphadamage, decay, r232Thdam', one(T), one(T))
    # Sm-147
    @. decay = exp(λ147Sm*(agesteps + dt_2)) - exp(λ147Sm*(agesteps - dt_2))
    mul!(alphadeposition, decay, r147SmHe', one(T), one(T))
    mul!(alphadamage, decay, r147Smdam', one(T), one(T))

    # Allocate additional variables that will be needed for Crank-Nicolson
    annealeddamage = similar(alphadamage)
    β = zeros(T, nrsteps)

    # Allocate arrays for diffusivities
    DL = zeros(T, length(tsteps))
    Dtrap = zeros(T, length(tsteps))

    # Allocate output matrix for all timesteps
    u = zeros(T, nrsteps, length(tsteps)+1)

    # Allocate variables for tridiagonal matrix and RHS
    dl = ones(T, nrsteps-1)    # Sub-diagonal row
    d = ones(T, nrsteps)       # Diagonal
    du = ones(T, nrsteps-1)    # Supra-diagonal row
    du2 = ones(T, nrsteps-2)   # sup-sup-diagonal row for pivoting

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = Tridiagonal(dl, d, du, du2)
    F = lu(A, allowsingular=true)

    # Vector for RHS of Crank-Nicolson equation with regular grid cells
    y = zeros(T, nrsteps)

    return ApatiteHe(
        T(age),
        T(age_sigma),
        T(offset),
        agesteps,
        tsteps,
        rsteps,
        redges,
        relvolumes,
        nrsteps,
        r238U,
        r235U,
        r232Th,
        r147Sm,
        alphadeposition,
        alphadamage,
        pr,
        annealeddamage,
        u,
        β,
        DL,
        Dtrap,
        A,
        F,
        y,
    )
end


"""
```julia
SphericalHe(T=Float64;
    age::Number = T(NaN),                   # [Ma] raw helium age
    age_sigma::Number = T(NaN),             # [Ma] raw helium age uncertainty (one-sigma)
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    stoppingpower::Number = T(1.189),       # [unitless] alpha stopping power relative to apatite
    r::Number,                              # [um] spherical radius 
    dr::Number = one(T),                    # [um] radial step size
    U238::Number,                           # [ppm] mineral U-238 concentration
    Th232::Number,                          # [ppm] mineral Th-232 concentration
    Sm147::Number = zero(T),                # [ppm] mineral Sm-147 concentration
    U238_matrix::Number = zero(T),          # [ppm] matrix U-238 concentration
    Th232_matrix::Number = zero(T),         # [ppm] matrix Th-232 concentration
    Sm147_matrix::Number = zero(T),         # [ppm] matrix Sm-147 concentration
    agesteps::AbstractVector | tsteps::AbstractVector, # Temporal discretization
)
```
Construct a `SphericalHe` chronometer representing a mineral with a raw 
helium age of `age` ± `age_sigma` [Ma], uniform diffusivity, a spherical
radius of `r` [μm], and uniform U, Th and Sm concentrations specified
by `U238`, `Th232`, and `Sm147` [PPMw]. (A present day U-235/U-238 
ratio of 1/137.818 is assumed)

Spatial discretization follows a halfwidth step of `dr` [μm], while temporal 
discretization  follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` 
must be sorted in increasing order.
"""
struct SphericalHe{T<:AbstractFloat} <: HeliumSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to other samples
    agesteps::FloatRange        # [Ma] age in Ma relative to the present
    tsteps::FloatRange          # [Ma] forward time since crystallization
    rsteps::FloatRange          # [um] radius bin centers
    redges::FloatRange          # [um] radius bin edges
    relvolumes::Vector{T}       # [unitless] fraction of volume in each radial step
    nrsteps::Int                # [n] number of radial steps, including both implicit points at each side
    r238U::Vector{T}            # [atoms/g] radial U-238 concentrations
    r235U::Vector{T}            # [atoms/g] radial U-235 concentrations
    r232Th::Vector{T}           # [atoms/g] radial Th-232 concentrations
    r147Sm::Vector{T}           # [atoms/g] radial Sm-147 concentrations
    alphadeposition::Matrix{T}  # [atoms/g] alpha deposition matrix
    u::Matrix{T}
    β::Vector{T}
    De::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{T, Tridiagonal{T, Vector{T}}, Vector{Int64}}
    y::Vector{T}
end
function SphericalHe(T::Type{<:AbstractFloat}=Float64;
        age::Number = T(NaN),
        age_sigma::Number = T(NaN),
        offset::Number = zero(T),
        stoppingpower::Number = T(1.189),
        r::Number, 
        dr::Number = one(T), 
        U238::Number, 
        Th232::Number, 
        Sm147::Number = zero(T), 
        U238_matrix::Number = zero(T), 
        Th232_matrix::Number = zero(T), 
        Sm147_matrix::Number = zero(T), 
        volumeweighting::Symbol = :spherical,
        agesteps = nothing,
        tsteps = nothing,
    )

    # Temporal discretization
    tsteps, agesteps = checkdiscretization(tsteps, agesteps)
    tsteps, agesteps = floatrange(tsteps), floatrange(agesteps)

    # Alpha stopping distances for each isotope in each decay chain, adjusted from those of apatite
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alpharadii238U = (13.54, 16.26, 15.84, 16.31, 20.09, 22.89, 33.39, 19.10,)./stoppingpower
    alpharadii235U = (14.48, 17.39, 22.5, 20.97, 26.89, 31.40, 26.18,)./stoppingpower
    alpharadii232Th = (12.60, 19.32, 21.08, 20.09, 27.53, 34.14,)./stoppingpower
    # Ketcham et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alpharadii147Sm = (5.93,)./stoppingpower

    # Crystal size and spatial discretization
    redges = floatrange(0 : dr : r)                 # Edges of each radius element
    rsteps = cntr(redges)                           # Centers of each radius element
    nrsteps = length(rsteps)+2                      # Number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    # Relative volume each concentric radial shell
    relvolumes = volumefraction(volumeweighting, redges, r)

    # Additional discretization outside of grain, for alpha injection
    outsideredges = floatrange(r : dr : r+maximum(alpharadii238U))
    outsidersteps = cntr(outsideredges)
    outsiderelvolumes = volumefraction(volumeweighting, outsideredges, r)

    # Observed radial HPE profiles at present day
    r238U = fill(T(U238), size(rsteps))         # [PPMw]
    r235U = fill(T(U238/137.818), size(rsteps)) # [PPMw]
    r232Th = fill(T(Th232), size(rsteps))       # [PPMw]
    r147Sm = fill(T(Sm147), size(rsteps))       # [PPMw]

    # Convert to atoms per gram
    r238U .*= 6.022E23 / 1E6 / 238
    r235U .*= 6.022E23 / 1E6 / 235
    r232Th .*= 6.022E23 / 1E6 / 232
    r147Sm .*= 6.022E23 / 1E6 / 147

    # Outside (bulk/matrix) HPE concentrations, in atoms per gram
    o238U = U238_matrix * 6.022E23 / 1E6 / 238
    o235U = U238_matrix/137.818 * 6.022E23 / 1E6 / 235
    o232Th = Th232_matrix * 6.022E23 / 1E6 / 232
    o147Sm = Sm147_matrix * 6.022E23 / 1E6 / 147

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance

    # Allocate intersection density vector
    dint = zeros(T, length(redges) - 1)

    # Radial alpha deposition from U-238
    r238UHe = zeros(T, size(r238U))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, relvolumes, r238U)
        for alpharadius in alpharadii238U
            sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,rsteps[ri])
            @. r238UHe += relvolumes[ri] * dint * r238U[ri]
        end
    end
    # Correct for alpha injection
    if o238U > 0
        @inbounds for ri in eachindex(outsidersteps, outsiderelvolumes)
            for alpharadius in alpharadii238U
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,outsidersteps[ri])
                @. r238UHe += outsiderelvolumes[ri] * dint * o238U
            end
        end
    end

    # Radial alpha deposition from U-235
    r235UHe = zeros(T, size(r235U))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, relvolumes, r235U)
        for alpharadius in alpharadii235U
            sphereintersectiondensity!(dint, redges,relvolumes,alpharadius,rsteps[ri])
            @. r235UHe += relvolumes[ri] * dint * r235U[ri]
        end
    end
    # Correct for alpha injection
    if o235U > 0
        @inbounds for ri in eachindex(outsidersteps, outsiderelvolumes)
            for alpharadius in alpharadii235U
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,outsidersteps[ri])
                @. r235UHe += outsiderelvolumes[ri] * dint * o235U
            end
        end
    end

    # Radial alpha deposition from Th-232
    r232ThHe = zeros(T, size(r232Th))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, relvolumes, r232Th)
        for alpharadius in alpharadii232Th
            sphereintersectiondensity!(dint, redges,relvolumes,alpharadius,rsteps[ri])
            @. r232ThHe += relvolumes[ri] * dint * r232Th[ri]
        end
    end
    # Correct for alpha injection
    if o232Th > 0
        @inbounds for ri in eachindex(outsidersteps, outsiderelvolumes)
            for alpharadius in alpharadii232Th
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,outsidersteps[ri])
                @. r232ThHe += outsiderelvolumes[ri] * dint * o232Th
            end
        end
    end

    # Radial alpha deposition from Sm-147
    r147SmHe = zeros(T, size(r147Sm))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, relvolumes, r147Sm)
        for alpharadius in alpharadii147Sm
            sphereintersectiondensity!(dint, redges,relvolumes,alpharadius,rsteps[ri])
            @. r147SmHe += relvolumes[ri] * dint * r147Sm[ri]
        end
    end
    # Correct for alpha injection
    if o147Sm > 0
        @inbounds for ri in eachindex(outsidersteps, outsiderelvolumes)
            for alpharadius in alpharadii147Sm
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                sphereintersectiondensity!(dint,redges,relvolumes,alpharadius,outsidersteps[ri])
                @. r147SmHe += outsiderelvolumes[ri] * dint * o147Sm
            end
        end
    end

    # Calculate corrected alpha deposition and recoil damage each time step for each radius
    dt_2 = step(tsteps)/2
    decay = zeros(T, length(tsteps))
    # Allocate deposition arrays
    alphadeposition = zeros(T, length(tsteps), nrsteps-2)

    # U-238
    @. decay = exp(λ238U*(agesteps + dt_2)) - exp(λ238U*(agesteps - dt_2))
    mul!(alphadeposition, decay, r238UHe', one(T), one(T))
    # U-235
    @. decay = exp(λ235U*(agesteps + dt_2)) - exp(λ235U*(agesteps - dt_2))
    mul!(alphadeposition, decay, r235UHe', one(T), one(T))
    # Th-232
    @. decay = exp(λ232Th*(agesteps + dt_2)) - exp(λ232Th*(agesteps - dt_2))
    mul!(alphadeposition, decay, r232ThHe', one(T), one(T))
    # Sm-147
    @. decay = exp(λ147Sm*(agesteps + dt_2)) - exp(λ147Sm*(agesteps - dt_2))
    mul!(alphadeposition, decay, r147SmHe', one(T), one(T))

    # Allocate additional variables that will be needed for Crank-Nicolson
    β = zeros(T, nrsteps)

    # Allocate arrays for diffusivities
    De = zeros(T, length(tsteps))

    # Allocate output matrix for all timesteps
    u = zeros(T, nrsteps, length(tsteps)+1)

    # Allocate variables for tridiagonal matrix and RHS
    dl = ones(T, nrsteps-1)    # Sub-diagonal row
    d = ones(T, nrsteps)       # Diagonal
    du = ones(T, nrsteps-1)    # Supra-diagonal row
    du2 = ones(T, nrsteps-2)   # sup-sup-diagonal row for pivoting

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = Tridiagonal(dl, d, du, du2)
    F = lu(A, allowsingular=true)

    # Vector for RHS of Crank-Nicolson equation with regular grid cells
    y = zeros(T, nrsteps)

    return SphericalHe(
        T(age),
        T(age_sigma),
        T(offset),
        agesteps,
        tsteps,
        rsteps,
        redges,
        relvolumes,
        nrsteps,
        r238U,
        r235U,
        r232Th,
        r147Sm,
        alphadeposition,
        u,
        β,
        De,
        A,
        F,
        y,
    )
end

"""
```julia
PlanarHe(T=Float64;
    age::Number = T(NaN),                   # [Ma] raw helium age
    age_sigma::Number = T(NaN),             # [Ma] raw helium age uncertainty (one-sigma)
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    stoppingpower::Number = T(1.189),       # [unitless] alpha stopping power relative to apatite
    r::Number,                              # [um] planar half-width
    dr::Number = one(T),                    # [um] radial step size
    U238::Number,                           # [ppm] mineral U-238 concentration
    Th232::Number,                          # [ppm] mineral Th-232 concentration
    Sm147::Number = zero(T),                # [ppm] mineral Sm-147 concentration
    U238_matrix::Number = zero(T),          # [ppm] matrix U-238 concentration
    Th232_matrix::Number = zero(T),         # [ppm] matrix Th-232 concentration
    Sm147_matrix::Number = zero(T),         # [ppm] matrix Sm-147 concentration
    agesteps::AbstractVector | tsteps::AbstractVector, # Temporal discretization
)
```
Construct an `PlanarHe` chronometer representing a mineral with a raw 
helium age of `age` ± `age_sigma` [Ma], uniform diffusivity, a planar
half-width of `r` [μm], and uniform U, Th and Sm concentrations specified
by `U238`, `Th232`, and `Sm147` [PPMw]. (A present day U-235/U-238 
ratio of 1/137.818 is assumed)

Spatial discretization follows a halfwidth step of `dr` [μm], while temporal 
discretization  follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` 
must be sorted in increasing order.
"""
struct PlanarHe{T<:AbstractFloat} <: HeliumSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to other samples
    agesteps::FloatRange        # [Ma] age in Ma relative to the present
    tsteps::FloatRange          # [Ma] forward time since crystallization
    rsteps::FloatRange          # [um] halfwidth bin centers
    redges::FloatRange          # [um] halfwidth bin edges
    nrsteps::Int                # [n] number of spatial steps, including both implicit points at each side
    r238U::Vector{T}            # [atoms/g] radial U-238 concentrations
    r235U::Vector{T}            # [atoms/g] radial U-235 concentrations
    r232Th::Vector{T}           # [atoms/g] radial Th-232 concentrations
    r147Sm::Vector{T}           # [atoms/g] radial Sm-147 concentrations
    alphadeposition::Matrix{T}  # [atoms/g] alpha deposition matrix
    u::Matrix{T}
    β::Vector{T}
    De::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{T, Tridiagonal{T, Vector{T}}, Vector{Int64}}
    y::Vector{T}
end
function PlanarHe(T::Type{<:AbstractFloat}=Float64;
        age::Number = T(NaN),
        age_sigma::Number = T(NaN),
        offset::Number = zero(T),
        stoppingpower::Number = T(1.189),
        r::Number, 
        dr::Number = one(T), 
        U238::Number, 
        Th232::Number, 
        Sm147::Number = zero(T), 
        U238_matrix::Number = zero(T), 
        Th232_matrix::Number = zero(T), 
        Sm147_matrix::Number = zero(T), 
        agesteps = nothing,
        tsteps = nothing,
    )

    # Temporal discretization
    tsteps, agesteps = checkdiscretization(tsteps, agesteps)
    tsteps, agesteps = floatrange(tsteps), floatrange(agesteps)

    # Alpha stopping distances for each isotope in each decay chain, adjusted from those of apatite
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alpharadii238U = (13.54, 16.26, 15.84, 16.31, 20.09, 22.89, 33.39, 19.10,)./stoppingpower
    alpharadii235U = (14.48, 17.39, 22.5, 20.97, 26.89, 31.40, 26.18,)./stoppingpower
    alpharadii232Th = (12.60, 19.32, 21.08, 20.09, 27.53, 34.14,)./stoppingpower
    # Ketcham et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alpharadii147Sm = (5.93,)./stoppingpower

    # Crystal size and spatial discretization
    redges = floatrange(0 : dr : r)                 # Edges of each halfwidth element
    rsteps = cntr(redges)                           # Centers of each halfwidth element
    nrsteps = length(rsteps)+2                      # Number of radial grid points -- note 2 implict points: one at negative halfwidth, one outside grain
    
    # Additional discretization outside of grain, for alpha injection
    outsideredges = floatrange(r : dr : r+maximum(alpharadii238U))
    outsidersteps = cntr(outsideredges)

    # Observed radial HPE profiles at present day
    r238U = fill(T(U238), size(rsteps))         # [PPMw]
    r235U = fill(T(U238/137.818), size(rsteps)) # [PPMw]
    r232Th = fill(T(Th232), size(rsteps))       # [PPMw]
    r147Sm = fill(T(Sm147), size(rsteps))       # [PPMw]

    # Convert to atoms per gram
    r238U .*= 6.022E23 / 1E6 / 238
    r235U .*= 6.022E23 / 1E6 / 235
    r232Th .*= 6.022E23 / 1E6 / 232
    r147Sm .*= 6.022E23 / 1E6 / 147

    # Outside (bulk/matrix) HPE concentrations, in atoms per gram
    o238U = U238_matrix * 6.022E23 / 1E6 / 238
    o235U = U238_matrix/137.818 * 6.022E23 / 1E6 / 235
    o232Th = Th232_matrix * 6.022E23 / 1E6 / 232
    o147Sm = Sm147_matrix * 6.022E23 / 1E6 / 147

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance

    # Allocate intersection density vector
    dint = zeros(T, length(redges) - 1)

    # Radial alpha deposition from U-238
    r238UHe = zeros(T, size(r238U))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, r238U)
        for alpharadius in alpharadii238U
            slabsphereintersectiondensity!(dint, redges,alpharadius,rsteps[ri])
            @. r238UHe += dint * r238U[ri]
        end
    end
    # Correct for alpha injection
    if o238U > 0
        @inbounds for ri in eachindex(outsidersteps)
            for alpharadius in alpharadii238U
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                slabsphereintersectiondensity!(dint, redges,alpharadius,outsidersteps[ri])
                @. r238UHe += dint * o238U
            end
        end
    end

    # Radial alpha deposition from U-235
    r235UHe = zeros(T, size(r235U))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, r235U)
        for alpharadius in alpharadii235U
            slabsphereintersectiondensity!(dint, redges,alpharadius,rsteps[ri])
            @. r235UHe += dint * r235U[ri]
        end
    end
    # Correct for alpha injection
    if o235U > 0
        @inbounds for ri in eachindex(outsidersteps)
            for alpharadius in alpharadii235U
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                slabsphereintersectiondensity!(dint, redges,alpharadius,outsidersteps[ri])
                @. r235UHe += dint * o235U
            end
        end
    end

    # Radial alpha deposition from Th-232
    r232ThHe = zeros(T, size(r232Th))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, r232Th)
        for alpharadius in alpharadii232Th
            slabsphereintersectiondensity!(dint, redges,alpharadius,rsteps[ri])
            @. r232ThHe += dint * r232Th[ri]
        end
    end
    # Correct for alpha injection
    if o232Th > 0
        @inbounds for ri in eachindex(outsidersteps)
            for alpharadius in alpharadii232Th
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                slabsphereintersectiondensity!(dint, redges,alpharadius,outsidersteps[ri])
                @. r232ThHe += dint * o232Th
            end
        end
    end

    # Radial alpha deposition from Sm-147
    r147SmHe = zeros(T, size(r147Sm))
    # Correct for alpha ejection
    @inbounds for ri in eachindex(rsteps, r147Sm)
        for alpharadius in alpharadii147Sm
            slabsphereintersectiondensity!(dint, redges,alpharadius,rsteps[ri])
            @. r147SmHe += dint * r147Sm[ri]
        end
    end
    # Correct for alpha injection
    if o147Sm > 0
        @inbounds for ri in eachindex(outsidersteps)
            for alpharadius in alpharadii147Sm
                (outsidersteps[ri] - first(outsidersteps)) > alpharadius && continue
                slabsphereintersectiondensity!(dint, redges,alpharadius,outsidersteps[ri])
                @. r147SmHe += dint * o147Sm
            end
        end
    end

    # Calculate corrected alpha deposition and recoil damage each time step for each halfwidth
    dt_2 = step(tsteps)/2
    decay = zeros(T, length(tsteps))
    # Allocate deposition arrays
    alphadeposition = zeros(T, length(tsteps), nrsteps-2)

    # U-238
    @. decay = exp(λ238U*(agesteps + dt_2)) - exp(λ238U*(agesteps - dt_2))
    mul!(alphadeposition, decay, r238UHe', one(T), one(T))
    # U-235
    @. decay = exp(λ235U*(agesteps + dt_2)) - exp(λ235U*(agesteps - dt_2))
    mul!(alphadeposition, decay, r235UHe', one(T), one(T))
    # Th-232
    @. decay = exp(λ232Th*(agesteps + dt_2)) - exp(λ232Th*(agesteps - dt_2))
    mul!(alphadeposition, decay, r232ThHe', one(T), one(T))
    # Sm-147
    @. decay = exp(λ147Sm*(agesteps + dt_2)) - exp(λ147Sm*(agesteps - dt_2))
    mul!(alphadeposition, decay, r147SmHe', one(T), one(T))

    # Allocate additional variables that will be needed for Crank-Nicolson
    β = zeros(T, nrsteps)

    # Allocate arrays for diffusivities
    De = zeros(T, length(tsteps))

    # Allocate output matrix for all timesteps
    u = zeros(T, nrsteps, length(tsteps)+1)

    # Allocate variables for tridiagonal matrix and RHS
    dl = ones(T, nrsteps-1)    # Sub-diagonal row
    d = ones(T, nrsteps)       # Diagonal
    du = ones(T, nrsteps-1)    # Supra-diagonal row
    du2 = ones(T, nrsteps-2)   # sup-sup-diagonal row for pivoting

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = Tridiagonal(dl, d, du, du2)
    F = lu(A, allowsingular=true)

    # Vector for RHS of Crank-Nicolson equation with regular grid cells
    y = zeros(T, nrsteps)

    return PlanarHe(
        T(age),
        T(age_sigma),
        T(offset),
        agesteps,
        tsteps,
        rsteps,
        redges,
        nrsteps,
        r238U,
        r235U,
        r232Th,
        r147Sm,
        alphadeposition,
        u,
        β,
        De,
        A,
        F,
        y,
    )
end

"""
```julia
SphericalAr(T=Float64;
    age::Number = T(NaN),                   # [Ma] Ar-40/Ar-39 age
    age_sigma::Number = T(NaN),             # [Ma] Ar-40/Ar-39 age uncertainty (one-sigma)
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    r::Number,                              # [um] equivalent spherical radius 
    dr::Number = one(T),                    # [um] radial step size
    K40::Number=16.34,                      # [ppm] mineral K-40 concentration
    agesteps::AbstractVector | tsteps::AbstractVector, # Temporal discretization
)
```
Construct an `SphericalAr` chronometer representing a mineral with a raw 
argon age of `age` ± `age_sigma` [Ma], a uniform diffusivity, a spherical
radius of `r` [μm], and uniform K-40 concentrations specified by `K40` [PPMw].

Spatial discretization follows a halfwidth step of `dr` [μm], while temporal 
discretization  follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` 
must be sorted in increasing order.
"""
struct SphericalAr{T<:AbstractFloat} <: ArgonSample{T}
    age::T                      # [Ma] Ar-40/Ar-39 age
    age_sigma::T                # [Ma] Ar-40/Ar-39 age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to other samples
    agesteps::FloatRange        # [Ma] age in Ma relative to the present
    tsteps::FloatRange          # [Ma] forward time since crystallization
    rsteps::FloatRange          # [um] radius bin centers
    redges::FloatRange          # [um] radius bin edges
    relvolumes::Vector{T}       # [unitless] fraction of volume in each radial step
    nrsteps::Int                # [n] number of radial steps, including both implicit points at each side
    r40K::Vector{T}             # [atoms/g] radial K-40 concentrations
    argondeposition::Matrix{T}  # [atoms/g] Ar-40 deposition matrix
    step_parent::Vector{T}
    step_daughter::Vector{T}
    u::Matrix{T}
    β::Vector{T}
    De::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{T, Tridiagonal{T, Vector{T}}, Vector{Int64}}
    y::Vector{T}
end
function SphericalAr(T::Type{<:AbstractFloat}=Float64;
        age::Number = T(NaN),
        age_sigma::Number = T(NaN),
        offset::Number = zero(T),
        r::Number, 
        dr::Number = one(T), 
        K40::Number = 16.34, 
        volumeweighting::Symbol = :spherical,
        agesteps = nothing,
        tsteps = nothing,
    )

    # Temporal discretization
    tsteps, agesteps = checkdiscretization(tsteps, agesteps)
    tsteps, agesteps = floatrange(tsteps), floatrange(agesteps)

    # Crystal size and spatial discretization
    rsteps = floatrange(0+dr/2 : dr : r-dr/2)
    redges = floatrange(     0 : dr : r     )   # Edges of each radius element
    nrsteps = length(rsteps)+2                  # number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    # Relative volume each concentric radial shell
    relvolumes = volumefraction(volumeweighting, redges, r)

    # Observed radial HPE profiles at present day
    r40K = fill(T(K40), size(rsteps))         # [PPMw]

    # Convert to atoms per gram
    r40K .*= 6.022E23 / 1E6 / 39.96399848
    # The proportion of that which will decay to Ar
    r40KAr = r40K .* BR40K

    # Calculate corrected argon deposition and recoil damage each time step for each radius
    dt_2 = step(tsteps)/2
    decay = zeros(T, length(tsteps))
    # Allocate deposition arrays
    argondeposition = zeros(T, length(tsteps), nrsteps-2)

    # K-40
    @. decay = exp(λ40K*(agesteps + dt_2)) - exp(λ40K*(agesteps - dt_2))
    mul!(argondeposition, decay, r40KAr', one(T), one(T))

    # Allocate additional variables that will be needed for Crank-Nicolson
    β = zeros(T, nrsteps)

    # Allocate arrays for diffusivities
    De = zeros(T, length(tsteps))

    # Allocate arrays to optionaly track parent and daughter concentrations
    step_parent = zeros(T, length(tsteps))
    step_daughter = zeros(T, length(tsteps))

    # Allocate output matrix for all timesteps
    u = zeros(T, nrsteps, length(tsteps)+1)

    # Allocate variables for tridiagonal matrix and RHS
    dl = ones(T, nrsteps-1)    # Sub-diagonal row
    d = ones(T, nrsteps)       # Diagonal
    du = ones(T, nrsteps-1)    # Supra-diagonal row
    du2 = ones(T, nrsteps-2)   # sup-sup-diagonal row for pivoting

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = Tridiagonal(dl, d, du, du2)
    F = lu(A, allowsingular=true)

    # Vector for RHS of Crank-Nicolson equation with regular grid cells
    y = zeros(T, nrsteps)

    return SphericalAr(
        T(age),
        T(age_sigma),
        T(offset),
        agesteps,
        tsteps,
        rsteps,
        redges,
        relvolumes,
        nrsteps,
        r40K,
        argondeposition,
        step_parent,
        step_daughter,
        u,
        β,
        De,
        A,
        F,
        y,
    )
end

"""
```julia
PlanarAr(T=Float64;
    age::Number = T(NaN),                   # [Ma] Ar-40/Ar-39 age
    age_sigma::Number = T(NaN),             # [Ma] Ar-40/Ar-39 age uncertainty (one-sigma)
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    r::Number,                              # [um] planar half-width
    dr::Number = one(T),                    # [um] radial step size
    K40::Number=16.34,                      # [ppm] mineral K-40 concentration
    agesteps::AbstractVector | tsteps::AbstractVector, # Temporal discretization
)
```
Construct an `PlanarAr` chronometer representing a mineral with a raw 
argon age of `age` ± `age_sigma` [Ma], a uniform diffusivity,
a radius of `r` [μm], and uniform K-40 concentrations specified by `K40` [PPM].

Spatial discretization follows a halfwidth step of `dr` [μm], while temporal 
discretization  follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` 
must be sorted in increasing order.
"""
struct PlanarAr{T<:AbstractFloat} <: ArgonSample{T}
    age::T                      # [Ma] Ar-40/Ar-39 age
    age_sigma::T                # [Ma] Ar-40/Ar-39 age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to other samples
    agesteps::FloatRange        # [Ma] age in Ma relative to the present
    tsteps::FloatRange          # [Ma] forward time since crystallization
    rsteps::FloatRange          # [um] halfwidth bin centers
    redges::FloatRange          # [um] halfwidth bin edges
    nrsteps::Int                # [n] number of spatial steps, including both implicit points at each side
    r40K::Vector{T}             # [atoms/g] radial K-40 concentrations
    argondeposition::Matrix{T}  # [atoms/g] Ar-40 deposition matrix
    step_parent::Vector{T}
    step_daughter::Vector{T}
    u::Matrix{T}
    β::Vector{T}
    De::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{T, Tridiagonal{T, Vector{T}}, Vector{Int64}}
    y::Vector{T}
end
function PlanarAr(T::Type{<:AbstractFloat}=Float64;
        age::Number = T(NaN),
        age_sigma::Number = T(NaN),
        offset::Number = zero(T),
        r::Number, 
        dr::Number = one(T), 
        K40::Number = 16.34, 
        agesteps = nothing,
        tsteps = nothing,
    )

    # Temporal discretization
    tsteps, agesteps = checkdiscretization(tsteps, agesteps)
    tsteps, agesteps = floatrange(tsteps), floatrange(agesteps)

    # Crystal size and spatial discretization
    rsteps = floatrange(0+dr/2 : dr : r-dr/2)
    redges = floatrange(     0 : dr : r     )   # Edges of each radius element
    nrsteps = length(rsteps)+2                  # number of radial grid points -- note 2 implict points: one at negative radius, one outside grain

    # Observed radial HPE profiles at present day
    r40K = fill(T(K40), size(rsteps))         # [PPMw]

    # Convert to atoms per gram
    r40K .*= 6.022E23 / 1E6 / 39.96399848
    # The proportion of that which will decay to Ar
    r40KAr = r40K .* BR40K

    # Calculate corrected argon deposition each time step for each radius
    dt_2 = step(tsteps)/2
    decay = zeros(T, length(tsteps))
    # Allocate deposition arrays
    argondeposition = zeros(T, length(tsteps), nrsteps-2)

    # K-40
    @. decay = exp(λ40K*(agesteps + dt_2)) - exp(λ40K*(agesteps - dt_2))
    mul!(argondeposition, decay, r40KAr', one(T), one(T))

    # Allocate additional variables that will be needed for Crank-Nicholson
    β = zeros(T, nrsteps)

    # Allocate arrays for diffusivities
    De = zeros(T, length(tsteps))

    # Allocate arrays to optionaly track parent and daughter concentrations
    step_parent = zeros(T, length(tsteps))
    step_daughter = zeros(T, length(tsteps))

    # Allocate output matrix for all timesteps
    u = zeros(T, nrsteps, length(tsteps)+1)

    # Allocate variables for tridiagonal matrix and RHS
    dl = ones(T, nrsteps-1)    # Sub-diagonal row
    d = ones(T, nrsteps)       # Diagonal
    du = ones(T, nrsteps-1)    # Supra-diagonal row
    du2 = ones(T, nrsteps-2)   # sup-sup-diagonal row for pivoting

    # Tridiagonal matrix for LHS of Crank-Nicholson equation with regular grid cells
    A = Tridiagonal(dl, d, du, du2)
    F = lu(A, allowsingular=true)

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = zeros(T, nrsteps)

    return PlanarAr(
        T(age),
        T(age_sigma),
        T(offset),
        agesteps,
        tsteps,
        rsteps,
        redges,
        nrsteps,
        r40K,
        argondeposition,
        step_parent,
        step_daughter,
        u,
        β,
        De,
        A,
        F,
        y,
    )
end

## --- Multiple domain diffusion chronometers!

    """
    ```julia
    MultipleDomain(T=Float64, C=PlanarAr;
        age::AbstractVector,                            # [Ma] measured Ar-40/Ar-39 ages at each degassing step
        age_sigma::AbstractVector,                      # [Ma] measured Ar-40/Ar-39 age uncertainties (one-sigma) at each degassing step
        fraction_experimental::AbstractVector,          # [unitless] cumulative fraction of total Ar-39 released each degassing step
        fraction_experimental_sigma::Number=T(0.01),    # [unitless] uncertainty in degassing fraction
        tsteps_experimental::AbstractVector,            # [s] time steps of experimental heating schedule
        Tsteps_experimental::AbstractVector,            # [C] temperature steps of experimental heating schedule
        fit::AbstractVector,                            # [Bool] Whether or not each degassing step should be used in inversion
        offset::Number = zero(T),                       # [C] temperature offset relative to other samples
        fuse::Bool = true,                              # [Bool] Treat the grain as having fused (released all remaining Ar)
        volume_fraction::AbstractVector,                # [unitless] fraction of total volume represented by each domain
        r::Number = 100,                                # [um] nominal model domain radius (spherical) or half-width (planar)
        dr::Number = one(T),                            # [um] nominal model domain radius step
        K40::Number = 16.34,                            # [ppm] mineral K-40 concentration
        agesteps::AbstractVector | tsteps::AbstractVector, # Temporal discretization
    )
    ```
    Construct a `MultipleDomain` diffusion chronometer given an observed argon
    release spectrum, degassing schedule, where each domain is represented by a 
    `PlanarAr` or `SphericalAr` chronometer.

    Domain diffusivity and volume parameters must be supplied as vectors
    `Ea` [kJ/mol], `lnD0a2` [log(1/s)], and `volume_fraction` [unitless]
    obtained by separately fitting the release spectrum (the former two
    as an `MDDiffusivity` object).

    See also: `MDDiffusivity`, `PlanarAr`, `SphericalAr`, `degas!`
    """
    struct MultipleDomain{T<:AbstractFloat, C<:Union{SphericalAr{T}, PlanarAr{T}}} <: AbsoluteChronometer{T}
        age::Vector{T}                      # [Ma] measured Ar-40/Ar-39 ages at each degassing step
        age_sigma::Vector{T}                # [Ma] measured Ar-40/Ar-39 age uncertainties (one-sigma) at each degassing step
        fraction_experimental::Vector{T}    # [unitless] cumulative fraction of total Ar-39 released each degassing step
        fraction_experimental_sigma::T      # [unitless] uncertainty in degassing fraction
        midpoint_experimental::Vector{T}    # [unitless] midpoint of fraction_experimental for each step
        tsteps_experimental::Vector{T}      # [s] time steps of experimental heating schedule
        Tsteps_experimental::Vector{T}      # [C] temperature steps of experimental heating schedule
        fit::BitVector                      # [Bool] Whether or not each step should be used in inversion
        offset::T                           # [C] temperature offset relative to other samples
        fuse::Bool                          # [Bool] Treat the grain as having fused (released all remaining Ar)
        domains::Vector{C}                  # Vector of chronometer obects for each domain
        volume_fraction::Vector{T}          # [unitless] fraction of total volume represented by each domain
        model_age::Vector{T}                # [Ma] calculated age at each model degassing step
        model_parent::Vector{T}             # [atoms/g equivalent] parent tracer degassed
        model_daughter::Vector{T}           # [atoms/g] daughter degassed
        model_fraction::Vector{T}           # [unitless] cumulative fraction of parent degasssed
        tsteps_degassing::FloatRange        # [s] time steps of model heating schedule
        Tsteps_degassing::Vector{T}         # [C] temperature steps of model heating schedule
        agesteps::FloatRange                # [Ma] age in Ma relative to the present
        tsteps::FloatRange                  # [Ma] forward time since crystallization
    end
    function MultipleDomain(T=Float64, C=PlanarAr;
            age::AbstractVector,
            age_sigma::AbstractVector,
            fraction_experimental::AbstractVector,
            fraction_experimental_sigma::Number=T(0.01),
            tsteps_experimental::AbstractVector,
            Tsteps_experimental::AbstractVector,
            fit::AbstractVector,
            offset::Number = zero(T),
            fuse::Bool = true,
            volume_fraction::AbstractVector,
            r::Number = 100,
            dr::Number = one(T),
            K40::Number = 16.34, 
            agesteps=nothing,
            tsteps=nothing,
        )
        # Temporal discretization
        isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
        isnothing(tsteps) && (tsteps = (maximum(agesteps)+minimum(agesteps)) .- agesteps)
        isnothing(agesteps) && (agesteps = (maximum(tsteps)+minimum(tsteps)) .- tsteps)
        agesteps, tsteps = floatrange(agesteps), floatrange(tsteps)
        @assert issorted(tsteps)
        
        # Check input arrays are the right size and ordered properly
        @assert eachindex(age) == eachindex(age_sigma) == eachindex(fraction_experimental) == eachindex(tsteps_experimental) == eachindex(Tsteps_experimental)
        @assert issorted(tsteps_experimental, lt=<=) "Degassing time steps must be in strictly increasing order"
        @assert all(x->0<=x<=1, fraction_experimental) "All \"fraction degassed\" values must be between 0 and 1"

        # Calculate midpoints of `fraction_experimental`
        midpoint_experimental = @. T(fraction_experimental + [0; fraction_experimental[1:end-1]])/2

        # Interpolate degassing t-T steps to the same resolution as the forward model
        tsteps_degassing = floatrange(first(tsteps_experimental), last(tsteps_experimental), length=length(agesteps))
        Tsteps_degassing = linterp1(tsteps_experimental, T.(Tsteps_experimental), tsteps_degassing) 
        model_age = zeros(T, length(tsteps_degassing))
        model_parent = zeros(T, length(tsteps_degassing))
        model_daughter = zeros(T, length(tsteps_degassing))
        model_fraction = zeros(T, length(tsteps_degassing))

        # Ensure volume fraction sums to one
        volume_fraction ./= nansum(volume_fraction) 
        
        # Allocate domains
        domains = [C(T; age=nanmean(age), age_sigma=nanstd(age), offset, r, dr, K40, agesteps) for i in eachindex(volume_fraction)]
        return MultipleDomain{T,C{T}}(
            T.(age),
            T.(age_sigma),
            T.(fraction_experimental),
            T(fraction_experimental_sigma),
            midpoint_experimental,
            tsteps_experimental,
            Tsteps_experimental,
            Bool.(fit),
            T(offset),
            fuse,
            domains,
            T.(volume_fraction),
            model_age,
            model_parent,
            model_daughter,
            model_fraction,
            tsteps_degassing,
            Tsteps_degassing,
            floatrange(agesteps),
            floatrange(tsteps),
        )
    end


## --- Check validity and consistency of temporal discretization

function checkdiscretization(tsteps, agesteps)
    isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
    isnothing(tsteps) && (tsteps = (maximum(agesteps)+minimum(agesteps)) .- agesteps)
    isnothing(agesteps) && (agesteps = (maximum(tsteps)+minimum(tsteps)) .- tsteps)
    @assert issorted(tsteps) "`tsteps` must be sorted in increasing order"
    @assert tsteps ≈ (maximum(agesteps)+minimum(agesteps)) .- agesteps "`tsteps` and `agesteps must represent the same chronology"
    return tsteps, agesteps
end

## --- Utility functions to get values, uncertainties, etc. from any chronometers

value(x::AbsoluteChronometer{T}) where {T} = x.age::T
stdev(x::AbsoluteChronometer{T}) where {T} = x.age_sigma::T
value(x::MultipleDomain{T}) where {T} = nanmean(x.age, @.(x.fit/x.age_sigma^2))::T
stdev(x::MultipleDomain{T}) where {T} = nanstd(x.age, @.(x.fit/x.age_sigma^2))::T
value(x::FissionTrackLength{T}) where {T} = x.length::T
stdev(x::FissionTrackLength{T}) where {T} = zero(T)
value(x::ApatiteTrackLengthOriented{T}) where {T} = x.lcmod::T

function val(x::Chronometer)
    @warn "Thermochron.val has been deprecated in favor of Thermochron.value"
    value(x)
end
function err(x::Chronometer)
    @warn "Thermochron.err has been deprecated in favor of Thermochron.stdev"
    stdev(x)
end

temperatureoffset(x::Chronometer) = x.offset

eU(x::Chronometer{T}) where {T<:AbstractFloat} = T(NaN)
function eU(x::HeliumSample{T}) where {T<:AbstractFloat}
    # Convert from atoms/g to ppm
    eu = nanmean(x.r238U) / (6.022E23 / 1E6 / 238)
    eu += 0.238*nanmean(x.r232Th) / (6.022E23 / 1E6 / 232)
    eu += 0.012*nanmean(x.r147Sm) / (6.022E23 / 1E6 / 147)
    return T(eu)
end

## -- Utility functions related to age and age uncertinty of absolute chronometers

# Get age and age sigma from a vector of chronometers
function get_age(x::AbstractArray{<:Chronometer{T}}, ::Type{C}=AbsoluteChronometer{T}) where {T<:AbstractFloat, C<:AbsoluteChronometer}
    result = sizehint!(T[], length(x))
    for xᵢ in x
        if isa(xᵢ, C)
            push!(result, value(xᵢ))
        end
    end
    return result
end
function get_age_sigma(x::AbstractArray{<:Chronometer{T}}, ::Type{C}=AbsoluteChronometer{T}) where {T<:AbstractFloat, C<:AbsoluteChronometer}
    result = sizehint!(T[], length(x))
    for xᵢ in x
        if isa(xᵢ, C)
            push!(result, stdev(xᵢ))
        end
    end
    return result
end

## --- End of File
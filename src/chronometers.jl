## --- Fission track length types

"""
```julia
ApatiteTrackLengthOriented(T::Type{<:AbstractFloat}=Float64; 
    length::Number = NaN,                   # [um] fission track length
    angle::Number = NaN,                    # [degrees] track angle from the c-axis
    lcmod::Number = lcmod(length, angle),   # [um] model length of an equivalent c-axis parallel rack
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    l0::Number = NaN,                       # [um] Initial track length
    l0_sigma::Number = NaN,                 # [um] Initial track length unertainty
    dpar::Number = NaN,                     # [um] diameter parallel to track
    F::Number = NaN,                        # [APFU] F concentration, in atoms per formula unit
    Cl::Number = NaN,                       # [APFU] Cl concentration, in atoms per formula unit
    OH::Number = NaN,                       # [APFU] OH concentration, in atoms per formula unit
    rmr0::Number = NaN,                     # [unitless] annealing parameter
    agesteps::AbstracVector | tsteps::AbstracVector, # Temporal discretization
)
```
Construct an `ApatiteTrackLengthOriented` chronometer representing a single apatite fission 
track `length` um long, oriented at `angle` degrees to the c-axis, with a relative annealing  
resistance specified by `rmr0`, optionally at a constant temperature offset (relative 
to other samples) of `offset` [C].

If not provided directly, `rmr0` will be calculated, in order of preference:
1. from `F`, `Cl`, and `OH` together, via the `apatite_rmr0model` function
2. from `Cl` alone, via the `apatite_rmr0fromcl` function
3. from `dpar`, via the `apatite_rmr0fromdpar` functions
4. using a default fallback value of 0.83, if none of the above are provided.

If not provided directly, `l0` and `l0_sigma` will be estimated using the 
`apatite_l0fromdpar` function, where dpar in turn is estimated using the
`apatite_dparfromrmr0` function if not provided directly.

Temporal discretization follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` must be sorted 
in increasing order.
"""
struct ApatiteTrackLengthOriented{T<:AbstractFloat, V<:AbstractVector{T}} <: FissionTrackLength{T}
    length::T               # [um] track length
    angle::T                # [degrees] track angle from the c-axis
    lcmod::T                # [um] model length of an equivalent c-axis parallel rack
    offset::T               # [C] temperature offset relative to other samples
    l0::T                   # [um] Initial track length
    l0_sigma::T             # [um] Initial track length unertainty
    rmr0::T                 # [unitless] relative resistance to annealing (0=most, 1=least)
    r::Vector{T}            # [unitless] reduced track lengths for each timestep
    pr::Vector{T}           # [unitless] reduced track densities for each timestep
    calc::Vector{T}         # [um] last calculated mean and standard deviation
    agesteps::V             # [Ma] age in Ma relative to the present
    tsteps::V               # [Ma] forward time since crystallization
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
        agesteps = nothing, 
        tsteps = nothing, 
        r = zeros(T, size(agesteps)),
        pr = zeros(T, size(agesteps)),
        calc = zeros(T, 2),
    )
    # Temporal discretization
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)
    l0, l0_sigma, rmr0 = parseaftparams(;l0, l0_sigma, dpar, F, Cl, OH, rmr0, oriented=true)
    @assert eachindex(r) == eachindex(pr) == eachindex(agesteps) == eachindex(tsteps)
    h = hash((ApatiteTrackLengthOriented, offset, agesteps, l0, l0_sigma, rmr0))
    ApatiteTrackLengthOriented(
        T(length),
        T(angle),
        T(lcmod),
        T(offset),
        T(l0),
        T(l0_sigma),
        T(rmr0),
        r,
        pr,
        calc,
        agesteps,
        tsteps,
    )
end


"""
```julia
ApatiteTrackLength(T::Type{<:AbstractFloat}=Float64; 
    length::Number = NaN,                   # [um] fission track length
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    l0::Number = NaN,                       # [um] Initial track length
    l0_sigma::Number = NaN,                 # [um] Initial track length unertainty
    dpar::Number = NaN,                     # [um] diameter parallel to track
    F::Number = NaN,                        # [APFU] F concentration, in atoms per formula unit
    Cl::Number = NaN,                       # [APFU] Cl concentration, in atoms per formula unit
    OH::Number = NaN,                       # [APFU] OH concentration, in atoms per formula unit
    rmr0::Number = NaN,                     # [unitless] annealing parameter
    agesteps::AbstracVector | tsteps::AbstracVector, # Temporal discretization
)
```
Construct an `ApatiteTrackLengthOriented` chronometer representing a single apatite 
fission track `length` um long with a relative annealing resistance specified by `rmr0`, 
optionally at a constant temperature offset (relative to other samples) of `offset` [C].

If not provided directly, `rmr0` will be calculated, in order of preference:
1. from `F`, `Cl`, and `OH` together, via the `apatite_rmr0model` function
2. from `Cl` alone, via the `apatite_rmr0fromcl` function
3. from `dpar`, via the `apatite_rmr0fromdpar` functions
4. using a default fallback value of 0.83, if none of the above are provided.

If not provided directly, `l0` and `l0_sigma` will be estimated using the 
`apatite_l0fromdpar` function, where dpar in turn is estimated using the
`apatite_dparfromrmr0` function if not provided directly.

Temporal discretization follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` must be sorted 
in increasing order.
"""
struct ApatiteTrackLength{T<:AbstractFloat, V<:AbstractVector{T}} <: FissionTrackLength{T}
    length::T               # [um] track length
    offset::T               # [C] temperature offset relative to other samples
    l0::T                   # [um] Initial track length
    l0_sigma::T             # [um] Initial track length unertainty
    rmr0::T                 # [unitless] relative resistance to annealing (0=most, 1=least)
    r::Vector{T}            # [unitless] reduced track lengths for each timestep
    pr::Vector{T}           # [unitless] reduced track densities for each timestep
    calc::Vector{T}         # [um] last calculated mean and standard deviation
    agesteps::V             # [Ma] age in Ma relative to the present
    tsteps::V               # [Ma] forward time since crystallization
    hash::UInt64            # [unitless] unique ID for linking track lengths with the same offset, agesteps, l0, etc.
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
        agesteps = nothing, 
        tsteps = nothing, 
        r = zeros(T, size(agesteps)),
        pr = zeros(T, size(agesteps)),
        calc = zeros(T, 2),
    )
    # Temporal discretization
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)
    # Multikinetic fission track parameters
    l0, l0_sigma, rmr0 = parseaftparams(;l0, l0_sigma, dpar, F, Cl, OH, rmr0, oriented=false)
    @assert eachindex(r) == eachindex(pr) == eachindex(agesteps) == eachindex(tsteps)
    h = hash((ApatiteTrackLength, offset, agesteps, l0, l0_sigma, rmr0))
    ApatiteTrackLength(
        T(length),
        T(offset),
        T(l0),
        T(l0_sigma),
        T(rmr0),
        r,
        pr,
        calc,
        agesteps,
        tsteps,
        h,
    )
end

"""
```julia
ZirconTrackLength(T::Type{<:AbstractFloat}=Float64; 
    length::Number = NaN,                   # [um] fission track length
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    l0::Number = 11.17,                     # [um] Initial track length
    l0_sigma::Number = 0.051,               # [um] Initial track length unertainty    
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
struct ZirconTrackLength{T<:AbstractFloat, V<:AbstractVector{T}} <: FissionTrackLength{T}
    length::T               # [um] track length
    offset::T               # [C] temperature offset relative to other samples
    l0::T                   # [um] initial track length
    l0_sigma::T             # [um] initial track length uncertainty
    r::Vector{T}            # [unitless] reduced track lengths for each timestep
    pr::Vector{T}           # [unitless] reduced track densities for each timestep
    calc::Vector{T}         # [um] last calculated mean and standard deviation
    agesteps::V             # [Ma] age in Ma relative to the present
    tsteps::V               # [Ma] forward time since crystallization
    hash::UInt64            # [unitless] unique ID for linking track lengths with the same offset, agesteps, l0, etc.
end
function ZirconTrackLength(T::Type{<:AbstractFloat}=Float64; 
        length::Number = NaN, 
        offset::Number = zero(T),
        l0::Number = 11.17,
        l0_sigma::Number = 0.051,
        agesteps = nothing, 
        tsteps = nothing,
        r = zeros(T, size(agesteps)),
        pr = zeros(T, size(agesteps)),
        calc = zeros(T, 2),
    )
    # Temporal discretization
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)
    # Initial track length and uncertainty
    if isnan(l0) 
        l0 = 11.17
    end
    if isnan(l0_sigma)
        l0_sigma = 0.051
    end
    @assert eachindex(r) == eachindex(pr) == eachindex(agesteps) == eachindex(tsteps)
    h = hash((ZirconTrackLength, offset, agesteps, l0, l0_sigma))
    ZirconTrackLength(
        T(length),
        T(offset),
        T(l0),
        T(l0_sigma),
        r,
        pr,
        calc,
        agesteps,
        tsteps,
        h,
    )
end


"""
```julia
MonaziteTrackLength(T::Type{<:AbstractFloat} = Float64; 
    length::Number = NaN,                   # [um] fission track length
    offset::Number = zero(T),               # [C] temperature offset relative to other samples
    l0::Number = 10.60,                     # [um] Initial track length
    l0_sigma::Number = 0.19,                # [um] Initial track length unertainty    
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
struct MonaziteTrackLength{T<:AbstractFloat, V<:AbstractVector{T}} <: FissionTrackLength{T}
    length::T               # [um] track length
    offset::T               # [C] temperature offset relative to other samples
    l0::T                   # [um] initial track length
    l0_sigma::T             # [um] initial track length uncertainty
    r::Vector{T}            # [unitless] reduced track lengths for each timestep
    pr::Vector{T}           # [unitless] reduced track densities for each timestep
    calc::Vector{T}         # [um] last calculated mean and standard deviation
    agesteps::V             # [Ma] age in Ma relative to the present
    tsteps::V               # [Ma] forward time since crystallization
    hash::UInt64            # [unitless] unique ID for linking track lengths with the same offset, agesteps, l0, etc.
end
function MonaziteTrackLength(T::Type{<:AbstractFloat}=Float64; 
        length::Number = NaN,
        offset::Number = zero(T),
        l0::Number = 10.60,
        l0_sigma::Number = 0.19,
        agesteps = nothing, 
        tsteps = nothing, 
        r = zeros(T, size(agesteps)),
        pr = zeros(T, size(agesteps)),
        calc = zeros(T, 2),
    )
    # Temporal discretization
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)
    # Initial track length and uncertainty
    if isnan(l0) 
        l0 = 10.60
    end
    if isnan(l0_sigma)
        l0_sigma = 0.19
    end
    @assert eachindex(r) == eachindex(pr) == eachindex(agesteps) == eachindex(tsteps)
    h = hash((MonaziteTrackLength, offset, agesteps, l0, l0_sigma))
    MonaziteTrackLength(
        T(length),
        T(offset),
        T(l0),
        T(l0_sigma),
        r,
        pr,
        calc,
        agesteps,
        tsteps,
        h,
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
struct ZirconFT{T<:AbstractFloat, V<:AbstractVector{T}} <: FissionTrackSample{T}
    age::T                  # [Ma] fission track age
    age_sigma::T            # [Ma] fission track age uncertainty (one-sigma)
    offset::T               # [C] temperature offset relative to other samples
    agesteps::V             # [Ma] age in Ma relative to the present
    tsteps::V               # [Ma] forward time since crystallization
end
function ZirconFT(T::Type{<:AbstractFloat}=Float64; 
        age::Number = NaN,              # [Ma] fission track age
        age_sigma::Number = NaN,        # [Ma] fission track age uncertainty
        offset::Number = zero(T),       # [C] temperature offset relative to other samples
        agesteps = nothing, 
        tsteps = nothing, 
    )
    # Temporal discretization
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)
    ZirconFT(
        T(age),
        T(age_sigma),
        T(offset),
        agesteps,
        tsteps,
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
struct MonaziteFT{T<:AbstractFloat, V<:AbstractVector{T}} <: FissionTrackSample{T}
    age::T                  # [Ma] fission track age
    age_sigma::T            # [Ma] fission track age uncertainty (one-sigma)
    offset::T               # [C] temperature offset relative to other samples
    agesteps::V             # [Ma] age in Ma relative to the present
    tsteps::V               # [Ma] forward time since crystallization
end
function MonaziteFT(T::Type{<:AbstractFloat}=Float64; 
        age::Number = NaN,              # [Ma] fission track age
        age_sigma::Number = NaN,        # [Ma] fission track age uncertainty
        offset::Number = zero(T),       # [C] temperature offset relative to other samples
        agesteps=nothing,
        tsteps=nothing, 
    )
    # Temporal discretization
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)
    MonaziteFT(
        T(age),
        T(age_sigma),
        T(offset),
        agesteps,
        tsteps,
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
1. from `F`, `Cl`, and `OH` together, via the `apatite_rmr0model` function
2. from `Cl` alone, via the `apatite_rmr0fromcl` function
3. from `dpar`, via the `apatite_rmr0fromdpar` functions
4. using a default fallback value of 0.83, if none of the above are provided.

Temporal discretization follows the age steps specified by `agesteps` (age before present)
and/or `tsteps` (forward time since crystallization), in Ma, where `tsteps` must be sorted 
in increasing order.
"""
struct ApatiteFT{T<:AbstractFloat, V<:AbstractVector{T}} <: FissionTrackSample{T}
    age::T                  # [Ma] fission track age
    age_sigma::T            # [Ma] fission track age uncertainty (one-sigma)
    offset::T               # [C] temperature offset relative to other samples
    rmr0::T                 # [unitless] relative resistance to annealing (0=most, 1=least)
    agesteps::V             # [Ma] age in Ma relative to the present
    tsteps::V               # [Ma] forward time since crystallization
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
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)
    # Multikinetic fission track parameters
    if isnan(rmr0)
        s = F + Cl + OH
        rmr0 = if !isnan(s)
            apatite_rmr0model(F/s*2, Cl/s*2, OH/s*2)
        elseif !isnan(Cl)
            apatite_rmr0fromcl(Cl)
        elseif !isnan(dpar)
            apatite_rmr0fromdpar(dpar)
        else
            0.83
        end
    end
    ApatiteFT(
        T(age),
        T(age_sigma),
        T(offset),
        T(rmr0),
        agesteps,
        tsteps,
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
    grainsize_matrix::Number = one(T),      # [mm] average grain size of matrix rock
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
struct ZirconHe{T<:AbstractFloat, V<:AbstractVector{T}} <: HeliumSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to other samples
    agesteps::V                 # [Ma] age in Ma relative to the present
    tsteps::V                   # [Ma] forward time since crystallization
    rsteps::FloatRange          # [um] radius bin centers
    redges::FloatRange          # [um] radius bin edges
    relvolumes::Vector{T}       # [unitless] fraction of volume in each radial step
    nrsteps::Int                # [n] number of radial steps, including both implicit points at each side
    r238U::Vector{T}            # [atoms/g] radial U-238 concentrations
    r235U::Vector{T}            # [atoms/g] radial U-235 concentrations
    r232Th::Vector{T}           # [atoms/g] radial Th-232 concentrations
    r147Sm::Vector{T}           # [atoms/g] radial Sm-147 concentrations
    bulkgrainsize::T            # [mm] average grain size of the whole-rock matrix
    bulkdeposition::Vector{T}  # [atoms/g] alpha (helium) production outside grain
    deposition::Matrix{T}  # [atoms/g] alpha (helium) deposition matrix within grain
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
    step_tracer::Vector{T}      # [atoms/g] buffer for degassed He-3 when modelling experimental heating schedule
    step_daughter::Vector{T}    # [atoms/g] buffer for degassed He-4 when modelling experimental heating schedule
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
        grainsize_matrix::Number = one(T),
        volumeweighting::Symbol=:cylindrical,
        agesteps = nothing,
        tsteps = nothing,
    )

    # Temporal discretization
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)

    # Crystal size and spatial discretization
    redges = floatrange(0 : dr : r)                 # Edges of each radius element
    rsteps = cntr(redges)                           # Centers of each radius element
    nrsteps = length(rsteps)+2                      # Number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    # Relative volume each concentric radial shell
    relvolumes = volumefraction(volumeweighting, redges, r)

    # Observed radial HPE profiles at present day, in atoms per gram
    r238U = fill(T(U238 * 6.022e17 / 238), size(rsteps))          # [atoms/g], converted from PPMw
    r235U = fill(T(U238/137.818 * 6.022e17 / 235), size(rsteps))  # [atoms/g], converted from PPMw
    r232Th = fill(T(Th232 * 6.022e17 / 232), size(rsteps))        # [atoms/g], converted from PPMw
    r147Sm = fill(T(Sm147 * 6.022e17 / 147), size(rsteps))        # [atoms/g], converted from PPMw

    # Outside (bulk/matrix) HPE concentrations, in atoms per gram
    o238U = U238_matrix * 6.022e17 / 238                          # [atoms/g], converted from PPMw
    o235U = U238_matrix/137.818 * 6.022e17 / 235                  # [atoms/g], converted from PPMw
    o232Th = Th232_matrix * 6.022e17 / 232                        # [atoms/g], converted from PPMw
    o147Sm = Sm147_matrix * 6.022e17 / 147                        # [atoms/g], converted from PPMw

    # Zircon alpha stopping distances for each isotope in each decay chain, from
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alpharadii238U = (11.78, 14.09, 13.73, 14.13, 17.32, 16.69, 28.56, 16.48,)
    alpharadii235U = (12.58, 15.04, 19.36, 18.06, 23.07, 26.87, 22.47,)
    alpharadii232Th = (10.99, 16.67, 18.16, 17.32, 23.61, 29.19,)
    # Ketcham et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alpharadii147Sm = (4.76,)
    # Additional discretization outside of grain, for alpha injection
    redgesₒ = floatrange(r : dr : r+maximum(alpharadii238U))
    relvolumesₒ = volumefraction(volumeweighting, redgesₒ, r)

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance [parent per-decay concentration equivalents]
    r238UHe = alphacorrectionspherical(alpharadii238U, r238U, redges, relvolumes, o238U, redgesₒ, relvolumesₒ)
    r235UHe = alphacorrectionspherical(alpharadii235U, r235U, redges, relvolumes, o235U, redgesₒ, relvolumesₒ)
    r232ThHe = alphacorrectionspherical(alpharadii232Th, r232Th, redges, relvolumes, o232Th, redgesₒ, relvolumesₒ)
    r147SmHe = alphacorrectionspherical(alpharadii147Sm, r147Sm, redges, relvolumes, o147Sm, redgesₒ, relvolumesₒ)

    # Alpha decay recoil damage [parent per-decay concentration equivalents]
    r238Udam = 8*r238U # No smoothing of alpha damage, 8 alphas per 238U
    r235Udam = 7*r235U # No smoothing of alpha damage, 7 alphas per 235U
    r232Thdam = 6*r232Th # No smoothing of alpha damage, 6 alphas per 232 Th
    r147Smdam = 1*r147Sm # No smoothing of alpha damage, 1 alpha per 147 Sm

    # Calculate corrected alpha deposition and recoil damage each time step for each radius
    decay = zeros(T, length(tsteps))
    # Allocate deposition and damage arrays
    bulkgrainsize = T(grainsize_matrix)
    bulkdeposition = zeros(T, length(tsteps))
    deposition = zeros(T, length(tsteps), nrsteps-2)
    alphadamage = zeros(T, length(tsteps), nrsteps-2)
    pr = zeros(T, length(tsteps), length(tsteps))

    # Calculate bin edges given agesteps, assuming step boundaries are halfway between steps
    agebinedges = [agesteps[1]-step_at(agesteps, 1)/2; cntr(agesteps); agesteps[end]+step_at(agesteps, lastindex(agesteps))/2]
    leftedges, rightedges = agebinedges[1:end-1], agebinedges[2:end]

    # U-238
    @. decay = exp(λ238U*leftedges) - exp(λ238U*rightedges)
    mul!(bulkdeposition, decay, o238U, 8, one(T))
    mul!(deposition, decay, r238UHe', one(T), one(T))
    mul!(alphadamage, decay, r238Udam', one(T), one(T))
    # U-235
    @. decay = exp(λ235U*leftedges) - exp(λ235U*rightedges)
    mul!(bulkdeposition, decay, o235U, 7, one(T))
    mul!(deposition, decay, r235UHe', one(T), one(T))
    mul!(alphadamage, decay, r235Udam', one(T), one(T))
    # Th-232
    @. decay = exp(λ232Th*leftedges) - exp(λ232Th*rightedges)
    mul!(bulkdeposition, decay, o232Th, 6, one(T))
    mul!(deposition, decay, r232ThHe', one(T), one(T))
    mul!(alphadamage, decay, r232Thdam', one(T), one(T))
    # Sm-147
    @. decay = exp(λ147Sm*leftedges) - exp(λ147Sm*rightedges)
    mul!(bulkdeposition, decay, o147Sm, one(T), one(T))
    mul!(deposition, decay, r147SmHe', one(T), one(T))
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

    # Allocate arrays to optionaly track tracer and daughter concentrations during degassing
    step_tracer = zeros(T, length(tsteps))
    step_daughter = zeros(T, length(tsteps))

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
        bulkgrainsize,
        bulkdeposition,
        deposition,
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
        step_tracer,
        step_daughter,
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
    grainsize_matrix::Number = one(T),      # [mm] average grain size of matrix rock
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
struct ApatiteHe{T<:AbstractFloat, V<:AbstractVector{T}} <: HeliumSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to other samples
    agesteps::V                 # [Ma] age in Ma relative to the present
    tsteps::V                   # [Ma] forward time since crystallization
    rsteps::FloatRange          # [um] radius bin centers
    redges::FloatRange          # [um] radius bin edges
    relvolumes::Vector{T}       # [unitless] fraction of volume in each radial step
    nrsteps::Int                # [n] number of radial steps, including both implicit points at each side
    r238U::Vector{T}            # [atoms/g] radial U-238 concentrations
    r235U::Vector{T}            # [atoms/g] radial U-235 concentrations
    r232Th::Vector{T}           # [atoms/g] radial Th-232 concentrations
    r147Sm::Vector{T}           # [atoms/g] radial Sm-147 concentrations
    bulkgrainsize::T            # [mm] average grain size of the whole-rock matrix
    bulkdeposition::Vector{T}  # [atoms/g] alpha (helium) production outside grain
    deposition::Matrix{T}  # [atoms/g] alpha (helium) deposition matrix within grain
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
    step_tracer::Vector{T}      # [atoms/g] buffer for degassed He-3 when modelling experimental heating schedule
    step_daughter::Vector{T}    # [atoms/g] buffer for degassed He-4 when modelling experimental heating schedule
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
        grainsize_matrix::Number = one(T),
        volumeweighting::Symbol = :cylindrical,
        agesteps = nothing,
        tsteps = nothing,
    )

    # Temporal discretization
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)

    # Crystal size and spatial discretization
    redges = floatrange(0 : dr : r)                 # Edges of each radius element
    rsteps = cntr(redges)                           # Centers of each radius element
    nrsteps = length(rsteps)+2                      # Number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    # Relative volume each concentric radial shell
    relvolumes = volumefraction(volumeweighting, redges, r)

    # Observed radial HPE profiles at present day, in atoms per gram
    r238U = fill(T(U238 * 6.022e17 / 238), size(rsteps))          # [atoms/g], converted from PPMw
    r235U = fill(T(U238/137.818 * 6.022e17 / 235), size(rsteps))  # [atoms/g], converted from PPMw
    r232Th = fill(T(Th232 * 6.022e17 / 232), size(rsteps))        # [atoms/g], converted from PPMw
    r147Sm = fill(T(Sm147 * 6.022e17 / 147), size(rsteps))        # [atoms/g], converted from PPMw

    # Outside (bulk/matrix) HPE concentrations, in atoms per gram
    o238U = U238_matrix * 6.022e17 / 238                          # [atoms/g], converted from PPMw
    o235U = U238_matrix/137.818 * 6.022e17 / 235                  # [atoms/g], converted from PPMw
    o232Th = Th232_matrix * 6.022e17 / 232                        # [atoms/g], converted from PPMw
    o147Sm = Sm147_matrix * 6.022e17 / 147                        # [atoms/g], converted from PPMw

    # Apatite alpha stopping distances for each isotope in each decay chain, from
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alpharadii238U = (13.54, 16.26, 15.84, 16.31, 20.09, 22.89, 33.39, 19.10,)
    alpharadii235U = (14.48, 17.39, 22.5, 20.97, 26.89, 31.40, 26.18,)
    alpharadii232Th = (12.60, 19.32, 21.08, 20.09, 27.53, 34.14,)
    # Ketcham et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alpharadii147Sm = (5.93,)
    # Additional discretization outside of grain, for alpha injection
    redgesₒ = floatrange(r : dr : r+maximum(alpharadii238U))
    relvolumesₒ = volumefraction(volumeweighting, redgesₒ, r)

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance [parent per-decay concentration equivalents]
    r238UHe = alphacorrectionspherical(alpharadii238U, r238U, redges, relvolumes, o238U, redgesₒ, relvolumesₒ)
    r235UHe = alphacorrectionspherical(alpharadii235U, r235U, redges, relvolumes, o235U, redgesₒ, relvolumesₒ)
    r232ThHe = alphacorrectionspherical(alpharadii232Th, r232Th, redges, relvolumes, o232Th, redgesₒ, relvolumesₒ)
    r147SmHe = alphacorrectionspherical(alpharadii147Sm, r147Sm, redges, relvolumes, o147Sm, redgesₒ, relvolumesₒ)

    # Alpha decay recoil damage [parent per-decay concentration equivalents]
    r238Udam = 8*r238U # No smoothing of alpha damage, 8 alphas per 238U
    r235Udam = 7*r235U # No smoothing of alpha damage, 7 alphas per 235U
    r232Thdam = 6*r232Th # No smoothing of alpha damage, 6 alphas per 232 Th
    r147Smdam = 1*r147Sm # No smoothing of alpha damage, 1 alpha per 147 Sm

    # Calculate corrected alpha deposition and recoil damage each time step for each radius
    decay = zeros(T, length(tsteps))
    # Allocate deposition and damage arrays
    bulkgrainsize = T(grainsize_matrix)
    bulkdeposition = zeros(T, length(tsteps))
    deposition = zeros(T, length(tsteps), nrsteps-2)
    alphadamage = zeros(T, length(tsteps), nrsteps-2)
    pr = zeros(T, length(tsteps), length(tsteps))

    # Calculate bin edges given agesteps, assuming step boundaries are halfway between steps
    agebinedges = [agesteps[1]-step_at(agesteps, 1)/2; cntr(agesteps); agesteps[end]+step_at(agesteps, lastindex(agesteps))/2]
    leftedges, rightedges = agebinedges[1:end-1], agebinedges[2:end]

    # U-238
    @. decay = exp(λ238U*leftedges) - exp(λ238U*rightedges)
    mul!(bulkdeposition, decay, o238U, 8, one(T))
    mul!(deposition, decay, r238UHe', one(T), one(T))
    mul!(alphadamage, decay, r238Udam', one(T), one(T))
    # U-235
    @. decay = exp(λ235U*leftedges) - exp(λ235U*rightedges)
    mul!(bulkdeposition, decay, o235U, 7, one(T))
    mul!(deposition, decay, r235UHe', one(T), one(T))
    mul!(alphadamage, decay, r235Udam', one(T), one(T))
    # Th-232
    @. decay = exp(λ232Th*leftedges) - exp(λ232Th*rightedges)
    mul!(bulkdeposition, decay, o232Th, 6, one(T))
    mul!(deposition, decay, r232ThHe', one(T), one(T))
    mul!(alphadamage, decay, r232Thdam', one(T), one(T))
    # Sm-147
    @. decay = exp(λ147Sm*leftedges) - exp(λ147Sm*rightedges)
    mul!(bulkdeposition, decay, o147Sm, one(T), one(T))
    mul!(deposition, decay, r147SmHe', one(T), one(T))
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

    # Allocate arrays to optionaly track tracer and daughter concentrations during degassing
    step_tracer = zeros(T, length(tsteps))
    step_daughter = zeros(T, length(tsteps))

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
        bulkgrainsize,
        bulkdeposition,
        deposition,
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
        step_tracer,
        step_daughter,
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
    grainsize_matrix::Number = one(T),      # [mm] average grain size of matrix rock
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
struct SphericalHe{T<:AbstractFloat, V<:AbstractVector{T}} <: HeliumSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to other samples
    agesteps::V                 # [Ma] age in Ma relative to the present
    tsteps::V                   # [Ma] forward time since crystallization
    rsteps::FloatRange          # [um] radius bin centers
    redges::FloatRange          # [um] radius bin edges
    relvolumes::Vector{T}       # [unitless] fraction of volume in each radial step
    nrsteps::Int                # [n] number of radial steps, including both implicit points at each side
    r238U::Vector{T}            # [atoms/g] radial U-238 concentrations
    r235U::Vector{T}            # [atoms/g] radial U-235 concentrations
    r232Th::Vector{T}           # [atoms/g] radial Th-232 concentrations
    r147Sm::Vector{T}           # [atoms/g] radial Sm-147 concentrations
    bulkgrainsize::T            # [mm] average grain size of the whole-rock matrix
    bulkdeposition::Vector{T}   # [atoms/g] alpha (helium) production outside grain
    deposition::Matrix{T}       # [atoms/g] alpha (helium) deposition matrix within grain
    u::Matrix{T}
    β::Vector{T}
    De::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{T, Tridiagonal{T, Vector{T}}, Vector{Int64}}
    y::Vector{T}
    step_tracer::Vector{T}      # [atoms/g] buffer for degassed He-3 when modelling experimental heating schedule
    step_daughter::Vector{T}    # [atoms/g] buffer for degassed He-4 when modelling experimental heating schedule
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
        grainsize_matrix::Number = one(T),
        volumeweighting::Symbol = :spherical,
        agesteps = nothing,
        tsteps = nothing,
    )

    # Temporal discretization
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)

    # Crystal size and spatial discretization
    redges = floatrange(0 : dr : r)                 # Edges of each radius element
    rsteps = cntr(redges)                           # Centers of each radius element
    nrsteps = length(rsteps)+2                      # Number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    # Relative volume each concentric radial shell
    relvolumes = volumefraction(volumeweighting, redges, r)

    # Observed radial HPE profiles at present day
    r238U = fill(T(U238 * 6.022e17 / 238), size(rsteps))          # [atoms/g], converted from PPMw
    r235U = fill(T(U238/137.818 * 6.022e17 / 235), size(rsteps))  # [atoms/g], converted from PPMw
    r232Th = fill(T(Th232 * 6.022e17 / 232), size(rsteps))        # [atoms/g], converted from PPMw
    r147Sm = fill(T(Sm147 * 6.022e17 / 147), size(rsteps))        # [atoms/g], converted from PPMw

    # Outside (bulk/matrix) HPE concentrations, in atoms per gram
    o238U = U238_matrix * 6.022e17 / 238                          # [atoms/g], converted from PPMw
    o235U = U238_matrix/137.818 * 6.022e17 / 235                  # [atoms/g], converted from PPMw
    o232Th = Th232_matrix * 6.022e17 / 232                        # [atoms/g], converted from PPMw
    o147Sm = Sm147_matrix * 6.022e17 / 147                        # [atoms/g], converted from PPMw

    # Alpha stopping distances for each isotope in each decay chain, adjusted from those of apatite
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alpharadii238U = (13.54, 16.26, 15.84, 16.31, 20.09, 22.89, 33.39, 19.10,)./stoppingpower
    alpharadii235U = (14.48, 17.39, 22.5, 20.97, 26.89, 31.40, 26.18,)./stoppingpower
    alpharadii232Th = (12.60, 19.32, 21.08, 20.09, 27.53, 34.14,)./stoppingpower
    # Ketcham et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alpharadii147Sm = (5.93,)./stoppingpower
    # Additional discretization outside of grain, for alpha injection
    redgesₒ = floatrange(r : dr : r+maximum(alpharadii238U))
    relvolumesₒ = volumefraction(volumeweighting, redgesₒ, r)

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance [parent per-decay concentration equivalents]
    r238UHe = alphacorrectionspherical(alpharadii238U, r238U, redges, relvolumes, o238U, redgesₒ, relvolumesₒ)
    r235UHe = alphacorrectionspherical(alpharadii235U, r235U, redges, relvolumes, o235U, redgesₒ, relvolumesₒ)
    r232ThHe = alphacorrectionspherical(alpharadii232Th, r232Th, redges, relvolumes, o232Th, redgesₒ, relvolumesₒ)
    r147SmHe = alphacorrectionspherical(alpharadii147Sm, r147Sm, redges, relvolumes, o147Sm, redgesₒ, relvolumesₒ)

    # Calculate corrected alpha deposition and recoil damage each time step for each radius
    decay = zeros(T, length(tsteps))
    # Allocate deposition arrays
    bulkgrainsize = T(grainsize_matrix)
    bulkdeposition = zeros(T, length(tsteps))
    deposition = zeros(T, length(tsteps), nrsteps-2)

    # Calculate bin edges given agesteps, assuming step boundaries are halfway between steps
    agebinedges = [agesteps[1]-step_at(agesteps, 1)/2; cntr(agesteps); agesteps[end]+step_at(agesteps, lastindex(agesteps))/2]
    leftedges, rightedges = agebinedges[1:end-1], agebinedges[2:end]

    # U-238
    @. decay = exp(λ238U*leftedges) - exp(λ238U*rightedges)
    mul!(bulkdeposition, decay, o238U, 8, one(T))
    mul!(deposition, decay, r238UHe', one(T), one(T))
    # U-235
    @. decay = exp(λ235U*leftedges) - exp(λ235U*rightedges)
    mul!(bulkdeposition, decay, o235U, 7, one(T))
    mul!(deposition, decay, r235UHe', one(T), one(T))
    # Th-232
    @. decay = exp(λ232Th*leftedges) - exp(λ232Th*rightedges)
    mul!(bulkdeposition, decay, o232Th, 6, one(T))
    mul!(deposition, decay, r232ThHe', one(T), one(T))
    # Sm-147
    @. decay = exp(λ147Sm*leftedges) - exp(λ147Sm*rightedges)
    mul!(bulkdeposition, decay, o147Sm, one(T), one(T))
    mul!(deposition, decay, r147SmHe', one(T), one(T))

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

    # Allocate arrays to optionaly track tracer and daughter concentrations during degassing
    step_tracer = zeros(T, length(tsteps))
    step_daughter = zeros(T, length(tsteps))

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
        bulkgrainsize,
        bulkdeposition,
        deposition,
        u,
        β,
        De,
        A,
        F,
        y,
        step_tracer,
        step_daughter,
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
    grainsize_matrix::Number = one(T),      # [mm] average grain size of matrix rock
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
struct PlanarHe{T<:AbstractFloat, V<:AbstractVector{T}} <: HeliumSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to other samples
    agesteps::V                 # [Ma] age in Ma relative to the present
    tsteps::V                   # [Ma] forward time since crystallization
    rsteps::FloatRange          # [um] halfwidth bin centers
    redges::FloatRange          # [um] halfwidth bin edges
    nrsteps::Int                # [n] number of spatial steps, including both implicit points at each side
    r238U::Vector{T}            # [atoms/g] radial U-238 concentrations
    r235U::Vector{T}            # [atoms/g] radial U-235 concentrations
    r232Th::Vector{T}           # [atoms/g] radial Th-232 concentrations
    r147Sm::Vector{T}           # [atoms/g] radial Sm-147 concentrations
    bulkgrainsize::T            # [mm] average grain size of the whole-rock matrix
    bulkdeposition::Vector{T}   # [atoms/g] alpha (helium) production outside grain
    deposition::Matrix{T}       # [atoms/g] alpha (helium) deposition matrix within grain
    u::Matrix{T}
    β::Vector{T}
    De::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{T, Tridiagonal{T, Vector{T}}, Vector{Int64}}
    y::Vector{T}
    step_tracer::Vector{T}      # [atoms/g] buffer for degassed He-3 when modelling experimental heating schedule
    step_daughter::Vector{T}    # [atoms/g] buffer for degassed He-4 when modelling experimental heating schedule
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
        grainsize_matrix::Number = one(T),
        agesteps = nothing,
        tsteps = nothing,
    )

    # Temporal discretization
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)

    # Crystal size and spatial discretization
    redges = floatrange(0 : dr : r)                 # Edges of each halfwidth element
    rsteps = cntr(redges)                           # Centers of each halfwidth element
    nrsteps = length(rsteps)+2                      # Number of radial grid points -- note 2 implict points: one at negative halfwidth, one outside grain

    # Observed radial HPE profiles at present day
    r238U = fill(T(U238 * 6.022e17 / 238), size(rsteps))          # [atoms/g], converted from PPMw
    r235U = fill(T(U238/137.818 * 6.022e17 / 235), size(rsteps))  # [atoms/g], converted from PPMw
    r232Th = fill(T(Th232 * 6.022e17 / 232), size(rsteps))        # [atoms/g], converted from PPMw
    r147Sm = fill(T(Sm147 * 6.022e17 / 147), size(rsteps))        # [atoms/g], converted from PPMw

    # Outside (bulk/matrix) HPE concentrations, in atoms per gram
    o238U = U238_matrix * 6.022e17 / 238                          # [atoms/g], converted from PPMw
    o235U = U238_matrix/137.818 * 6.022e17 / 235                  # [atoms/g], converted from PPMw
    o232Th = Th232_matrix * 6.022e17 / 232                        # [atoms/g], converted from PPMw
    o147Sm = Sm147_matrix * 6.022e17 / 147                        # [atoms/g], converted from PPMw

    # Alpha stopping distances for each isotope in each decay chain, adjusted from those of apatite
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alpharadii238U = (13.54, 16.26, 15.84, 16.31, 20.09, 22.89, 33.39, 19.10,)./stoppingpower
    alpharadii235U = (14.48, 17.39, 22.5, 20.97, 26.89, 31.40, 26.18,)./stoppingpower
    alpharadii232Th = (12.60, 19.32, 21.08, 20.09, 27.53, 34.14,)./stoppingpower
    # Ketcham et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alpharadii147Sm = (5.93,)./stoppingpower
    # Additional discretization outside of grain, for alpha injection
    redgesₒ = floatrange(r : dr : r+maximum(alpharadii238U))

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance [parent per-decay concentration equivalents]
    r238UHe = alphacorrectionslab(alpharadii238U, r238U, redges, o238U, redgesₒ)
    r235UHe = alphacorrectionslab(alpharadii235U, r235U, redges, o235U, redgesₒ)
    r232ThHe = alphacorrectionslab(alpharadii232Th, r232Th, redges, o232Th, redgesₒ)
    r147SmHe = alphacorrectionslab(alpharadii147Sm, r147Sm, redges, o147Sm, redgesₒ)

    # Calculate corrected alpha deposition and recoil damage each time step for each halfwidth
    decay = zeros(T, length(tsteps))
    # Allocate deposition arrays
    bulkgrainsize = T(grainsize_matrix)
    bulkdeposition = zeros(T, length(tsteps))
    deposition = zeros(T, length(tsteps), nrsteps-2)

    # Calculate bin edges given agesteps, assuming step boundaries are halfway between steps
    agebinedges = [agesteps[1]-step_at(agesteps, 1)/2; cntr(agesteps); agesteps[end]+step_at(agesteps, lastindex(agesteps))/2]
    leftedges, rightedges = agebinedges[1:end-1], agebinedges[2:end]

    # U-238
    @. decay = exp(λ238U*leftedges) - exp(λ238U*rightedges)
    mul!(bulkdeposition, decay, o238U, 8, one(T))
    mul!(deposition, decay, r238UHe', one(T), one(T))
    # U-235
    @. decay = exp(λ235U*leftedges) - exp(λ235U*rightedges)
    mul!(bulkdeposition, decay, o235U, 7, one(T))
    mul!(deposition, decay, r235UHe', one(T), one(T))
    # Th-232
    @. decay = exp(λ232Th*leftedges) - exp(λ232Th*rightedges)
    mul!(bulkdeposition, decay, o232Th, 6, one(T))
    mul!(deposition, decay, r232ThHe', one(T), one(T))
    # Sm-147
    @. decay = exp(λ147Sm*leftedges) - exp(λ147Sm*rightedges)
    mul!(bulkdeposition, decay, o147Sm, one(T), one(T))
    mul!(deposition, decay, r147SmHe', one(T), one(T))

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

    # Allocate arrays to optionaly track tracer and daughter concentrations during degassing
    step_tracer = zeros(T, length(tsteps))
    step_daughter = zeros(T, length(tsteps))

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
        bulkgrainsize,
        bulkdeposition,
        deposition,
        u,
        β,
        De,
        A,
        F,
        y,
        step_tracer,
        step_daughter,
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
    K40::Number = 16.34,                    # [ppm] mineral K-40 concentration
    K40_matrix::Number = zero(T)            # [ppm] matrix K-40 concentration
    grainsize_matrix::Number = one(T),      # [mm] average grain size of matrix rock
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
struct SphericalAr{T<:AbstractFloat, V<:AbstractVector{T}} <: ArgonSample{T}
    age::T                      # [Ma] Ar-40/Ar-39 age
    age_sigma::T                # [Ma] Ar-40/Ar-39 age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to other samples
    agesteps::V                 # [Ma] age in Ma relative to the present
    tsteps::V                   # [Ma] forward time since crystallization
    rsteps::FloatRange          # [um] radius bin centers
    redges::FloatRange          # [um] radius bin edges
    relvolumes::Vector{T}       # [unitless] fraction of volume in each radial step
    nrsteps::Int                # [n] number of radial steps, including both implicit points at each side
    r40K::Vector{T}             # [atoms/g] radial K-40 concentrations
    bulkgrainsize::T            # [mm] average grain size of the whole-rock matrix
    bulkdeposition::Vector{T}   # [atoms/g], converted from PPMw Ar-40 production outside grain
    deposition::Matrix{T}       # [atoms/g], converted from PPMw Ar-40 deposition matrix
    u::Matrix{T}
    β::Vector{T}
    De::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{T, Tridiagonal{T, Vector{T}}, Vector{Int64}}
    y::Vector{T}
    step_tracer::Vector{T}      # [atoms/g] buffer for degassed Ar-39 when modelling experimental heating schedule
    step_daughter::Vector{T}    # [atoms/g] buffer for degassed Ar-40 when modelling experimental heating schedule
end
function SphericalAr(T::Type{<:AbstractFloat}=Float64;
        age::Number = T(NaN),
        age_sigma::Number = T(NaN),
        offset::Number = zero(T),
        r::Number, 
        dr::Number = one(T), 
        K40::Number = 16.34, 
        K40_matrix::Number = zero(T),
        grainsize_matrix::Number = one(T),
        volumeweighting::Symbol = :spherical,
        agesteps = nothing,
        tsteps = nothing,
    )

    # Temporal discretization
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)

    # Crystal size and spatial discretization
    rsteps = floatrange(0+dr/2 : dr : r-dr/2)
    redges = floatrange(     0 : dr : r     )   # Edges of each radius element
    nrsteps = length(rsteps)+2                  # number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    # Relative volume each concentric radial shell
    relvolumes = volumefraction(volumeweighting, redges, r)

    # Observed radial HPE profiles at present day
    r40K = fill(T(K40), size(rsteps))         # [PPMw]

    # Convert to atoms per gram
    r40K .*= 6.022e17 / 39.96399848

    # The proportion of that which will decay to Ar
    r40KAr = r40K .* BR40K

    # Outside (bulk/matrix) HPE concentrations, in atoms per gram
    o40K = K40_matrix * 6.022e17 / 39.96399848

    # Calculate corrected argon deposition and recoil damage each time step for each radius
    decay = zeros(T, length(tsteps))
    # Allocate deposition arrays
    bulkgrainsize = T(grainsize_matrix)
    bulkdeposition = zeros(T, length(tsteps))
    deposition = zeros(T, length(tsteps), nrsteps-2)

    # Calculate bin edges given agesteps, assuming step boundaries are halfway between steps
    agebinedges = [agesteps[1]-step_at(agesteps, 1)/2; cntr(agesteps); agesteps[end]+step_at(agesteps, lastindex(agesteps))/2]
    leftedges, rightedges = agebinedges[1:end-1], agebinedges[2:end]

    # K-40
    @. decay = exp(λ40K*leftedges) - exp(λ40K*rightedges)
    mul!(bulkdeposition, decay, o40K, BR40K, one(T))
    mul!(deposition, decay, r40KAr', one(T), one(T))

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

    # Allocate arrays to optionaly track tracer and daughter concentrations during degassing
    step_tracer = zeros(T, length(tsteps))
    step_daughter = zeros(T, length(tsteps))

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
        bulkgrainsize,
        bulkdeposition,
        deposition,
        u,
        β,
        De,
        A,
        F,
        y,
        step_tracer,
        step_daughter,
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
    K40::Number = 16.34,                    # [ppm] mineral K-40 concentration
    K40_matrix::Number = zero(T)            # [ppm] matrix K-40 concentration
    grainsize_matrix::Number = one(T),      # [mm] average grain size of matrix rock
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
struct PlanarAr{T<:AbstractFloat, V<:AbstractVector{T}} <: ArgonSample{T}
    age::T                      # [Ma] Ar-40/Ar-39 age
    age_sigma::T                # [Ma] Ar-40/Ar-39 age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to other samples
    agesteps::V                 # [Ma] age in Ma relative to the present
    tsteps::V                   # [Ma] forward time since crystallization
    rsteps::FloatRange          # [um] halfwidth bin centers
    redges::FloatRange          # [um] halfwidth bin edges
    nrsteps::Int                # [n] number of spatial steps, including both implicit points at each side
    r40K::Vector{T}             # [atoms/g] radial K-40 concentrations
    bulkgrainsize::T            # [mm] average grain size of the whole-rock matrix
    bulkdeposition::Vector{T} # [atoms/g], converted from PPMw Ar-40 production outside grain
    deposition::Matrix{T}  # [atoms/g], converted from PPMw Ar-40 deposition matrix
    u::Matrix{T}
    β::Vector{T}
    De::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{T, Tridiagonal{T, Vector{T}}, Vector{Int64}}
    y::Vector{T}
    step_tracer::Vector{T}      # [atoms/g] buffer for degassed Ar-39 when modelling experimental heating schedule
    step_daughter::Vector{T}    # [atoms/g] buffer for degassed Ar-40 when modelling experimental heating schedule
end
function PlanarAr(T::Type{<:AbstractFloat}=Float64;
        age::Number = T(NaN),
        age_sigma::Number = T(NaN),
        offset::Number = zero(T),
        r::Number, 
        dr::Number = one(T), 
        K40::Number = 16.34, 
        K40_matrix::Number = zero(T),
        grainsize_matrix::Number = one(T),
        agesteps = nothing,
        tsteps = nothing,
    )

    # Temporal discretization
    agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)

    # Crystal size and spatial discretization
    rsteps = floatrange(0+dr/2 : dr : r-dr/2)
    redges = floatrange(     0 : dr : r     )   # Edges of each radius element
    nrsteps = length(rsteps)+2                  # number of radial grid points -- note 2 implict points: one at negative radius, one outside grain

    # Observed radial HPE profiles at present day
    r40K = fill(T(K40), size(rsteps))         # [PPMw]

    # Convert to atoms per gram
    r40K .*= 6.022e17 / 39.96399848

    # The proportion of that which will decay to Ar
    r40KAr = r40K .* BR40K

    # Outside (bulk/matrix) HPE concentrations, in atoms per gram
    o40K = K40_matrix * 6.022e17 / 39.96399848

    # Calculate corrected argon deposition each time step for each radius
    decay = zeros(T, length(tsteps))
    # Allocate deposition arrays
    bulkgrainsize = T(grainsize_matrix)
    bulkdeposition = zeros(T, length(tsteps))
    deposition = zeros(T, length(tsteps), nrsteps-2)

    # Calculate bin edges given agesteps, assuming step boundaries are halfway between steps
    agebinedges = [agesteps[1]-step_at(agesteps, 1)/2; cntr(agesteps); agesteps[end]+step_at(agesteps, lastindex(agesteps))/2]
    leftedges, rightedges = agebinedges[1:end-1], agebinedges[2:end]

    # K-40
    @. decay = exp(λ40K*leftedges) - exp(λ40K*rightedges)
    mul!(bulkdeposition, decay, o40K, BR40K, one(T))
    mul!(deposition, decay, r40KAr', one(T), one(T))

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

    # Allocate arrays to optionaly track tracer and daughter concentrations during degassing
    step_tracer = zeros(T, length(tsteps))
    step_daughter = zeros(T, length(tsteps))

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
        bulkgrainsize,
        bulkdeposition,
        deposition,
        u,
        β,
        De,
        A,
        F,
        y,
        step_tracer,
        step_daughter,
    )
end

## --- Single-domain diffusion chronometers (e.g. 4/3 He)

    """
    ```julia
    SingleDomain(T=Float64, C=ApatiteHe;
        step_age::AbstractVector,                       # [Ma] measured Ar-40/Ar-39 ages at each degassing step
        step_age_sigma::AbstractVector,                 # [Ma] measured Ar-40/Ar-39 age uncertainties (one-sigma) at each degassing step
        fraction_experimental::AbstractVector,          # [unitless] cumulative fraction of total Ar-39 released each degassing step
        fraction_experimental_sigma=fill(T(0.005), size(fraction_experimental)),     # [unitless] uncertainty in degassing fraction
        tsteps_experimental::AbstractVector,            # [s] time steps of experimental heating schedule
        Tsteps_experimental::AbstractVector,            # [C] temperature steps of experimental heating schedule
        fit::AbstractVector,                            # [Bool] Whether or not each degassing step should be used in inversion
        offset::Number = zero(T),                       # [C] temperature offset relative to other samples
        fuse::Bool = true,                              # [Bool] Treat the grain as having fused (released all remaining Ar)
        agesteps::AbstractVector | tsteps::AbstractVector, # Temporal discretization
        kwargs...                                       # Any further keyword arguments are forwarded to the constructor for the domain `C`
    )
    ```
    Construct a `SingleDomain` diffusion chronometer given an observed 
    experimental release spectrum and degassing schedule, where the  
    single modelled domain is represented by a Chronometer object 
    (e.g. `ApatiteHe` for He-4/He-3).

    See also: `degas!`
    """
    struct SingleDomain{T<:AbstractFloat, C<:Union{HeliumSample{T}, ArgonSample{T}}} <: AbsoluteChronometer{T}
        step_age::Vector{T}                          # [Ma or unitless] measured ages (for Ar-40/Ar-39) or Rstep/Rbulk ratios (for He-4/He-3) at each degassing step
        step_age_sigma::Vector{T}                    # [Ma or unitless] measured age (or ratio) uncertainties (one-sigma) at each degassing step
        fraction_experimental::Vector{T}        # [unitless] cumulative fraction of total tracer (Ar-39 or He-3) released each degassing step
        fraction_experimental_sigma::Vector{T}  # [unitless] uncertainty in degassing fraction
        midpoint_experimental::Vector{T}        # [unitless] midpoint of fraction_experimental for each step
        tsteps_experimental::Vector{T}          # [s] time steps of experimental heating schedule
        Tsteps_experimental::Vector{T}          # [C] temperature steps of experimental heating schedule
        fit::BitVector                          # [Bool] Whether or not each step should be used in inversion
        offset::T                               # [C] temperature offset relative to other samples
        fuse::Bool                              # [Bool] Treat the grain as having fused (released all remaining Ar)
        domain::C                               # Chronometer obect for diffusion domain
        model_age::Vector{T}                    # [Ma or ratio] calculated age at each model degassing step
        model_fraction::Vector{T}               # [unitless] cumulative fraction of tracer He-3 degasssed
        tsteps_degassing::FloatRange            # [s] time steps of model heating schedule
        Tsteps_degassing::Vector{T}             # [C] temperature steps of model heating schedule
    end
    function SingleDomain(T=Float64, C=ApatiteHe;
            step_age::AbstractVector,
            step_age_sigma::AbstractVector,
            fraction_experimental::AbstractVector,
            fraction_experimental_sigma::AbstractVector=fill(T(0.005), size(fraction_experimental)),
            tsteps_experimental::AbstractVector,
            Tsteps_experimental::AbstractVector,
            fit::AbstractVector,
            offset::Number = zero(T),
            fuse::Bool = true,
            agesteps=nothing,
            tsteps=nothing,
            kwargs...
        )
        # Temporal discretization
        agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)
        
        # Check input arrays are the right size and ordered properly
        @assert eachindex(step_age) == eachindex(step_age_sigma) == eachindex(fraction_experimental) == eachindex(tsteps_experimental) == eachindex(Tsteps_experimental)
        @assert issorted(tsteps_experimental, lt=<=) "Degassing time steps must be in strictly increasing order"
        @assert all(x->0<=x<=1, fraction_experimental) "All \"fraction degassed\" values must be between 0 and 1"

        # Calculate midpoints of `fraction_experimental`
        midpoint_experimental = @. T(fraction_experimental + [0; fraction_experimental[1:end-1]])/2

        # Interpolate degassing t-T steps to the same resolution as the forward model
        tsteps_degassing = floatrange(first(tsteps_experimental), last(tsteps_experimental), length=length(agesteps))
        Tsteps_degassing = linterp1(tsteps_experimental, T.(Tsteps_experimental), tsteps_degassing) 
        model_age = zeros(T, length(tsteps_degassing))
        model_tracer = zeros(T, length(tsteps_degassing))
        model_daughter = zeros(T, length(tsteps_degassing))
        model_fraction = zeros(T, length(tsteps_degassing))
        
        # Allocate domain
        domain = C(T; offset, agesteps, tsteps, kwargs...)
        return SingleDomain(
            T.(step_age),
            T.(step_age_sigma),
            T.(fraction_experimental),
            T.(fraction_experimental_sigma),
            midpoint_experimental,
            tsteps_experimental,
            Tsteps_experimental,
            Bool.(fit),
            T(offset),
            fuse,
            domain,
            model_age,
            model_fraction,
            tsteps_degassing,
            Tsteps_degassing,
        )
    end

## --- Multiple domain diffusion chronometers!

    """
    ```julia
    MultipleDomain(T=Float64, C=PlanarAr;
        step_age::AbstractVector,                       # [Ma] measured Ar-40/Ar-39 ages at each degassing step
        step_age_sigma::AbstractVector,                 # [Ma] measured Ar-40/Ar-39 age uncertainties (one-sigma) at each degassing step
        fraction_experimental::AbstractVector,          # [unitless] cumulative fraction of total Ar-39 released each degassing step
        fraction_experimental_sigma=fill(T(0.005), size(fraction_experimental)),     # [unitless] uncertainty in degassing fraction
        tsteps_experimental::AbstractVector,            # [s] time steps of experimental heating schedule
        Tsteps_experimental::AbstractVector,            # [C] temperature steps of experimental heating schedule
        fit::AbstractVector,                            # [Bool] Whether or not each degassing step should be used in inversion
        offset::Number = zero(T),                       # [C] temperature offset relative to other samples
        fuse::Bool = true,                              # [Bool] Treat the grain as having fused (released all remaining Ar)
        volume_fraction::AbstractVector,                # [unitless] fraction of total volume represented by each domain
        r::Number = 100,                                # [um] nominal model domain radius (spherical) or half-width (planar)
        dr::Number = one(T),                            # [um] nominal model domain radius step
        agesteps::AbstractVector | tsteps::AbstractVector, # Temporal discretization
        kwargs...                                       # Any further keyword arguments are forwarded to the constructor for the domain `C`
    )
    ```
    Construct a `MultipleDomain` diffusion chronometer given an observed
    release spectrum and degassing schedule, where each domain is in turn 
    represented by a Chronometer object (e.g. `PlanarAr`, `SphericalAr`).

    Domain diffusivity and volume parameters must be supplied as vectors
    `Ea` [kJ/mol], `lnD0a2` [log(1/s)], and `volume_fraction` [unitless]
    obtained by separately fitting the release spectrum (the former two
    as an `MDDiffusivity` object).

    See also: `MDDiffusivity`, `PlanarAr`, `SphericalAr`, `degas!`
    """
    struct MultipleDomain{T<:AbstractFloat, C<:Union{HeliumSample{T}, ArgonSample{T}}} <: AbsoluteChronometer{T}
        step_age::Vector{T}                     # [Ma or unitless] measured ages (for Ar-40/Ar-39) or Rstep/Rbulk ratios (for He-4/He-3) at each degassing step
        step_age_sigma::Vector{T}               # [Ma or unitless] measured age (or ratio) uncertainties at each degassing step
        fraction_experimental::Vector{T}        # [unitless] cumulative fraction of total tracer (Ar-39 or He-3) released each degassing step
        fraction_experimental_sigma::Vector{T}  # [unitless] uncertainty in degassing fraction
        midpoint_experimental::Vector{T}        # [unitless] midpoint of fraction_experimental for each step
        tsteps_experimental::Vector{T}          # [s] time steps of experimental heating schedule
        Tsteps_experimental::Vector{T}          # [C] temperature steps of experimental heating schedule
        fit::BitVector                          # [Bool] Whether or not each step should be used in inversion
        offset::T                               # [C] temperature offset relative to other samples
        fuse::Bool                              # [Bool] Treat the grain as having fused (released all remaining Ar)
        domains::Vector{C}                      # Vector of chronometer obects for each domain
        volume_fraction::Vector{T}              # [unitless] fraction of total volume represented by each domain
        model_age::Vector{T}                    # [Ma] calculated age at each model degassing step
        model_tracer::Vector{T}                 # [atoms/g equivalent] parent tracer degassed
        model_daughter::Vector{T}               # [atoms/g], converted from PPMw daughter degassed
        model_fraction::Vector{T}               # [unitless] cumulative fraction of parent tracer degasssed
        tsteps_degassing::FloatRange            # [s] time steps of model heating schedule
        Tsteps_degassing::Vector{T}             # [C] temperature steps of model heating schedule
    end
    function MultipleDomain(T=Float64, C=PlanarAr;
            step_age::AbstractVector,
            step_age_sigma::AbstractVector,
            fraction_experimental::AbstractVector,
            fraction_experimental_sigma::AbstractVector=fill(T(0.005), size(fraction_experimental)),
            tsteps_experimental::AbstractVector,
            Tsteps_experimental::AbstractVector,
            fit::AbstractVector,
            offset::Number = zero(T),
            fuse::Bool = true,
            volume_fraction::AbstractVector,
            r::Number = 100,
            dr::Number = one(T),
            agesteps=nothing,
            tsteps=nothing,
            kwargs...
        )
        # Temporal discretization
        agesteps, tsteps = checktimediscretization(T, agesteps, tsteps)
        
        # Check input arrays are the right size and ordered properly
        @assert eachindex(step_age) == eachindex(step_age_sigma) == eachindex(fraction_experimental) == eachindex(fraction_experimental_sigma) == eachindex(tsteps_experimental) == eachindex(Tsteps_experimental)
        @assert issorted(tsteps_experimental, lt=<=) "Degassing time steps must be in strictly increasing order"
        @assert all(x->0<=x<=1, fraction_experimental) "All \"fraction degassed\" values must be between 0 and 1"

        # Calculate midpoints of `fraction_experimental`
        midpoint_experimental = @. T(fraction_experimental + [0; fraction_experimental[1:end-1]])/2
        @assert eachindex(midpoint_experimental) == eachindex(fraction_experimental)

        # Interpolate degassing t-T steps to the same resolution as the forward model
        tsteps_degassing = floatrange(first(tsteps_experimental), last(tsteps_experimental), length=length(agesteps))
        Tsteps_degassing = linterp1(tsteps_experimental, T.(Tsteps_experimental), tsteps_degassing) 
        model_age = zeros(T, length(tsteps_degassing))
        model_tracer = zeros(T, length(tsteps_degassing))
        model_daughter = zeros(T, length(tsteps_degassing))
        model_fraction = zeros(T, length(tsteps_degassing))

        # Ensure volume fraction sums to one
        if !isapprox(nansum(volume_fraction), 1, atol=0.01)
            @warn "volume fractions $volume_fraction do not sum to 1"
            volume_fraction ./= nansum(volume_fraction)
        end
        
        # Allocate domains
        bulk_age = nanmean(step_age, @.(fit./step_age_sigma^2))
        bulk_age_sigma = nanstd(step_age, @.(fit./step_age_sigma^2))
        domains = [C(T; age=bulk_age, age_sigma=bulk_age_sigma, offset, r, dr, agesteps, tsteps, kwargs...) for i in eachindex(volume_fraction)]
        return MultipleDomain(
            T.(step_age),
            T.(step_age_sigma),
            T.(fraction_experimental),
            T.(fraction_experimental_sigma),
            midpoint_experimental,
            tsteps_experimental,
            Tsteps_experimental,
            Bool.(fit),
            T(offset),
            fuse,
            domains,
            T.(volume_fraction),
            model_age,
            model_tracer,
            model_daughter,
            model_fraction,
            tsteps_degassing,
            Tsteps_degassing,
        )
    end

## --- Functions on Chronometers which require types to have been defined

# Implement eltype methods to deal with chronometers which are wrapper types
Base.eltype(x::Chronometer) = typeof(x)
Base.eltype(x::MultipleDomain{T,C}) where {T,C} = C
Base.eltype(x::SingleDomain{T,C}) where {T,C} = C

# Retrive the nominal value (age, length, etc) of any Chronometer
value(x::AbsoluteChronometer{T}) where {T} = x.age::T
value(x::Union{SingleDomain{T}, MultipleDomain{T}}) where {T} = nanmean(x.step_age, @.(x.fit/x.step_age_sigma^2))::T
value(x::FissionTrackLength{T}) where {T} = x.length::T
value(x::ApatiteTrackLengthOriented{T}) where {T} = x.lcmod::T
function val(x::Chronometer)
    @warn "Thermochron.val has been deprecated in favor of Thermochron.value"
    value(x)
end

# Retrive the nominal 1-sigma uncertainty (in age, length, etc.) of any Chronometer
stdev(x::AbsoluteChronometer{T}) where {T} = x.age_sigma::T
stdev(x::Union{SingleDomain{T}, MultipleDomain{T}}) where {T} = nanstd(x.step_age, @.(x.fit/x.step_age_sigma^2))::T
stdev(x::FissionTrackLength{T}) where {T} = zero(T)
function err(x::Chronometer)
    @warn "Thermochron.err has been deprecated in favor of Thermochron.stdev"
    stdev(x)
end

# Retrive the temperature offset of any Chronometer
temperatureoffset(x::Chronometer{T}) where {T} = x.offset::T

# Retrive the temporal discretization of any Chronometer
agediscretization(x::Chronometer{T}) where {T} = x.agesteps::AbstractVector{T}
agediscretization(x::MultipleDomain) = agediscretization(first(x.domains))
agediscretization(x::SingleDomain) = agediscretization(x.domain)
timediscretization(x::Chronometer{T}) where {T} = x.tsteps::AbstractVector{T}
timediscretization(x::MultipleDomain) = timediscretization(first(x.domains))
timediscretization(x::SingleDomain) = timediscretization(x.domain)

# Retrive the eU ("effective uranium") of any Chronometer
eU(x::Chronometer{T}) where {T<:AbstractFloat} = T(NaN)
function eU(x::HeliumSample{T}) where {T<:AbstractFloat}
    # Convert from atoms/g to ppm
    eu = nanmean(x.r238U) / (6.022e17 / 238)
    eu += 0.238*nanmean(x.r232Th) / (6.022e17 / 232)
    eu += 0.012*nanmean(x.r147Sm) / (6.022e17 / 147)
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
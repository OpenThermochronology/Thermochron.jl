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

eU(x::Chronometer{T}) where {T<:AbstractFloat} = T(NaN)
function eU(x::HeliumSample{T}) where {T<:AbstractFloat}
    # Convert from atoms/g to ppm
    eu = nanmean(x.r238U) / (6.022E23 / 1E6 / 238)
    eu += 0.238*nanmean(x.r232Th) / (6.022E23 / 1E6 / 232)
    eu += 0.012*nanmean(x.r147Sm) / (6.022E23 / 1E6 / 147)
    return T(eu)
end

## --- Fission track length types

struct ApatiteTrackLength{T<:AbstractFloat} <: FissionTrackLength{T}
    length::T               # [um] track length
    angle::T                # [degrees] track angle from the c-axis
    lcmod::T                # [um] model length of an equivalent c-axis parallel rack
    offset::T               # [C] temperature offset relative to the surface
    agesteps::FloatRange    # [Ma] age in Ma relative to the present
    tsteps::FloatRange      # [Ma] forward time since crystallization
    r::Vector{T}            # [unitless]
    pr::Vector{T}           # [unitless]
    ledges::FloatRange      # [um] Length distribution edges
    ldist::Vector{T}        # [um] Length log likelihood
    rmr0::T                 # [unitless] relative resistance to annealing (0=most, 1=least)
end
function ApatiteTrackLength(T::Type{<:AbstractFloat}=Float64; 
        length = T(NaN), 
        angle = T(NaN), 
        lcmod = lcmod(length, angle),
        offset::Number = zero(T),
        ledges = (0:1.0:20),
        dpar = T(NaN), 
        F = T(NaN), 
        Cl = T(NaN), 
        OH = T(NaN), 
        rmr0 = T(NaN),
        agesteps = nothing, 
        tsteps = nothing, 
    )
    # Temporal discretization
    isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
    isnothing(tsteps) && (tsteps=reverse(agesteps))
    isnothing(agesteps) && (agesteps=reverse(tsteps))
    agesteps, tsteps = floatrange(agesteps), floatrange(tsteps)
    @assert issorted(tsteps)
    # Multikinetic fission track parameters
    if isnan(rmr0)
        s = F + Cl + OH
        rmr0 = if !isnan(s)
            rmr0model(F/s*2, Cl/s*2, OH/s*2)
        elseif !isnan(dpar)
            rmr0fromdpar(dpar)
        else
            0.83
        end
    end
    r=zeros(T, size(agesteps))
    pr=zeros(T, size(agesteps))
    ldist=zeros(T, size(ledges).-1)
    ApatiteTrackLength(
        T(length),
        T(angle),
        T(lcmod),
        T(offset),
        floatrange(agesteps),
        floatrange(tsteps),
        r,
        pr,
        floatrange(ledges),
        ldist,
        T(rmr0),
    )
end

struct ZirconTrackLength{T<:AbstractFloat} <: FissionTrackLength{T}
    length::T               # [um] track length
    offset::T               # [C] temperature offset relative to the surface
    agesteps::FloatRange    # [Ma] age in Ma relative to the present
    tsteps::FloatRange      # [Ma] forward time since crystallization
    r::Vector{T}            # [unitless]
    pr::Vector{T}           # [unitless]
    ledges::FloatRange      # [um] Length distribution edges
    ldist::Vector{T}        # [um] Length log likelihood
end
function ZirconTrackLength(T::Type{<:AbstractFloat}=Float64; 
        length = T(NaN), 
        offset::Number = zero(T),
        ledges = (0:1.0:20),
        agesteps = nothing, 
        tsteps = nothing,
    )
    # Temporal discretization
    isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
    isnothing(tsteps) && (tsteps=reverse(agesteps))
    isnothing(agesteps) && (agesteps=reverse(tsteps))
    agesteps, tsteps = floatrange(agesteps), floatrange(tsteps)
    @assert issorted(tsteps)

    r=zeros(T, size(agesteps))
    pr=zeros(T, size(agesteps))
    ldist=zeros(T, size(ledges).-1)
    ZirconTrackLength(
        T(length),
        T(offset),
        floatrange(agesteps),
        floatrange(tsteps),
        r,
        pr,
        floatrange(ledges),
        ldist,
    )
end

struct MonaziteTrackLength{T<:AbstractFloat} <: FissionTrackLength{T}
    length::T               # [um] track length
    offset::T               # [C] temperature offset relative to the surface
    agesteps::FloatRange    # [Ma] age in Ma relative to the present
    tsteps::FloatRange      # [Ma] forward time since crystallization
    r::Vector{T}            # [unitless]
    pr::Vector{T}           # [unitless]
    ledges::FloatRange      # [um] Length distribution edges
    ldist::Vector{T}        # [um] Length log likelihood
end
function MonaziteTrackLength(T::Type{<:AbstractFloat}=Float64; 
        length = T(NaN), 
        offset::Number = zero(T),
        ledges = (0:1.0:20),
        agesteps = nothing, 
        tsteps = nothing, 
    )
    # Temporal discretization
    isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
    isnothing(tsteps) && (tsteps=reverse(agesteps))
    isnothing(agesteps) && (agesteps=reverse(tsteps))
    agesteps, tsteps = floatrange(agesteps), floatrange(tsteps)
    @assert issorted(tsteps)

    r=zeros(T, size(agesteps))
    pr=zeros(T, size(agesteps))
    ldist=zeros(T, size(ledges).-1)
    MonaziteTrackLength(
        T(length),
        T(offset),
        floatrange(agesteps),
        floatrange(tsteps),
        r,
        pr,
        floatrange(ledges),
        ldist,
    )
end

## --- Fission track age types

struct ZirconFT{T<:AbstractFloat} <: FissionTrackSample{T}
    age::T                  # [Ma] fission track age
    age_sigma::T            # [Ma] fission track age uncertainty (one-sigma)
    offset::T               # [C] temperature offset relative to the surface
    agesteps::FloatRange    # [Ma] age in Ma relative to the present
    tsteps::FloatRange      # [Ma] forward time since crystallization
end
function ZirconFT(T::Type{<:AbstractFloat}=Float64; 
        age = T(NaN), 
        age_sigma = T(NaN), 
        offset::Number = zero(T),
        agesteps = nothing, 
        tsteps = nothing, 
    )
    # Temporal discretization
    isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
    isnothing(tsteps) && (tsteps=reverse(agesteps))
    isnothing(agesteps) && (agesteps=reverse(tsteps))
    agesteps, tsteps = floatrange(agesteps), floatrange(tsteps)
    @assert issorted(tsteps)
    ZirconFT(
        T(age),
        T(age_sigma),
        T(offset),
        floatrange(agesteps),
        floatrange(tsteps),
    )
end

struct MonaziteFT{T<:AbstractFloat} <: FissionTrackSample{T}
    age::T                  # [Ma] fission track age
    age_sigma::T            # [Ma] fission track age uncertainty (one-sigma)
    offset::T               # [C] temperature offset relative to the surface
    agesteps::FloatRange    # [Ma] age in Ma relative to the present
    tsteps::FloatRange      # [Ma] forward time since crystallization
end
function MonaziteFT(T::Type{<:AbstractFloat}=Float64; 
        age = T(NaN), 
        age_sigma = T(NaN), 
        offset::Number = zero(T),
        agesteps, 
        tsteps=reverse(agesteps), 
    )
    MonaziteFT(
        T(age),
        T(age_sigma),
        T(offset),
        floatrange(agesteps),
        floatrange(tsteps),
    )
end

struct ApatiteFT{T<:AbstractFloat} <: FissionTrackSample{T}
    age::T                  # [Ma] fission track age
    age_sigma::T            # [Ma] fission track age uncertainty (one-sigma)
    offset::T               # [C] temperature offset relative to the surface
    agesteps::FloatRange    # [Ma] age in Ma relative to the present
    tsteps::FloatRange      # [Ma] forward time since crystallization
    rmr0::T                 # [unitless] relative resistance to annealing (0=most, 1=least)
end
function ApatiteFT(T::Type{<:AbstractFloat}=Float64; 
        age = T(NaN), 
        age_sigma = T(NaN),
        offset::Number = zero(T),
        agesteps, 
        tsteps=reverse(agesteps), 
        dpar = T(NaN),
        F = T(NaN),
        Cl = T(NaN),
        OH = T(NaN),
        rmr0 = T(NaN),
    )
    if isnan(rmr0)
        s = F + Cl + OH
        rmr0 = if !isnan(s)
            rmr0model(F/s*2, Cl/s*2, OH/s*2)
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
    age::Number = T(NaN),
    age_sigma::Number = T(NaN),
    T(offset),
    r::Number = one(T),
    dr::Number, 
    U238::Number,
    Th232::Number,
    Sm147::Number = zero(T),
    U238_matrix::Number = zero(T), 
    Th232_matrix::Number = zero(T), 
    Sm147_matrix::Number = zero(T), 
    volumeweighting::Symbol=:cylindrical,
    agesteps::AbstractRange,
)
```
Construct a `ZirconHe` chronometer representing a zircon with a raw 
helium age of `age` ± `age_sigma` [Ma], a  radius of `r` [μm], and uniform 
U, Th and Sm concentrations specified by `U238`, `Th232`, and `Sm147` [PPMw]. 
A present day U-235/U-238 ratio of 1/137.818 is assumed.

Spatial discretization follows a radius step of `dr` μm, and temporal
discretization follows the age steps specified by the `agesteps` range,
in Ma.
"""
struct ZirconHe{T<:AbstractFloat} <: HeliumSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to the surface
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
        r::Number = one(T),
        dr::Number, 
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
    isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
    isnothing(tsteps) && (tsteps=reverse(agesteps))
    isnothing(agesteps) && (agesteps=reverse(tsteps))
    agesteps, tsteps = floatrange(agesteps), floatrange(tsteps)
    @assert issorted(tsteps)

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
    volumeweighting::Symbol = :spherical,
    agesteps::AbstractRange,
)
```
Construct an `ApatiteHe` chronometer representing an apatite with a raw 
helium age of `age` ± `age_sigma` [Ma], a  radius of `r` [μm], and uniform 
U, Th and Sm concentrations specified by `U238`, `Th232`, and `Sm147` [PPMW]. 
A present day U-235/U-238 ratio of 1/137.818 is assumed.

Spatial discretization follows a radius step of `dr` [μm], and temporal
discretization follows the age steps specified by the `agesteps` range,
in Ma.
"""
struct ApatiteHe{T<:AbstractFloat} <: HeliumSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to the surface
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
    isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
    isnothing(tsteps) && (tsteps=reverse(agesteps))
    isnothing(agesteps) && (agesteps=reverse(tsteps))
    agesteps, tsteps = floatrange(agesteps), floatrange(tsteps)
    @assert issorted(tsteps)

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
    agesteps::AbstractRange,
)
```
Construct a `SphericalHe` chronometer representing a mineral with a raw 
helium age of `age` ± `age_sigma` [Ma], uniform diffusivity,
a radius of `r` [μm], and uniform U, Th and Sm concentrations specified
by `U238`, `Th232`, and `Sm147` [PPM]. (A present day U-235/U-238 
ratio of 1/137.818 is assumed)

Spatial discretization follows a radius step of `dr` [μm], and temporal
discretization follows the age steps specified by the `agesteps` range,
in Ma.
"""
struct SphericalHe{T<:AbstractFloat} <: HeliumSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to the surface
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
    isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
    isnothing(tsteps) && (tsteps=reverse(agesteps))
    isnothing(agesteps) && (agesteps=reverse(tsteps))
    agesteps, tsteps = floatrange(agesteps), floatrange(tsteps)
    @assert issorted(tsteps)

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
    agesteps::AbstractRange,
)
```
Construct an `PlanarHe` chronometer representing a mineral with a raw 
helium age of `age` ± `age_sigma` [Ma], uniform diffusivity, 
a halfwidth of `r` [μm], and uniform U, Th and Sm concentrations specified
by `U238`, `Th232`, and `Sm147` [PPM]. (A present day U-235/U-238 
ratio of 1/137.818 is assumed)

Spatial discretization follows a halfwidth step of `dr` [μm], and temporal
discretization follows the age steps specified by the `agesteps` range,
in Ma.
"""
struct PlanarHe{T<:AbstractFloat} <: HeliumSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to the surface
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
    isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
    isnothing(tsteps) && (tsteps=reverse(agesteps))
    isnothing(agesteps) && (agesteps=reverse(tsteps))
    agesteps, tsteps = floatrange(agesteps), floatrange(tsteps)
    @assert issorted(tsteps)

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
    age::Number = T(NaN),
    age_sigma::Number = T(NaN),
    offset::Number = zero(T),
    r::Number, 
    dr::Number = one(T), 
    K40::Number=16.34, 
    agesteps::AbstractRange,
)
```
Construct an `SphericalAr` chronometer representing a mineral with a raw 
argon age of `age` ± `age_sigma` [Ma], a uniform diffusivity,
a radius of `r` [μm], and uniform K-40 concentrations specified by `K40` [PPM].

Spatial discretization follows a radius step of `dr` [μm], and temporal
discretization follows the age steps specified by the `agesteps` range,
in Ma.
"""
struct SphericalAr{T<:AbstractFloat} <: ArgonSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to the surface
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
    isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
    isnothing(tsteps) && (tsteps=reverse(agesteps))
    isnothing(agesteps) && (agesteps=reverse(tsteps))
    agesteps, tsteps = floatrange(agesteps), floatrange(tsteps)
    @assert issorted(tsteps)

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
    age::Number = T(NaN),
    age_sigma::Number = T(NaN),
    offset::Number = zero(T),
    r::Number, 
    dr::Number = one(T), 
    K40::Number=16.34, 
    agesteps::AbstractRange,
)
```
Construct an `PlanarAr` chronometer representing a mineral with a raw 
argon age of `age` ± `age_sigma` [Ma], a uniform diffusivity,
a radius of `r` [μm], and uniform K-40 concentrations specified by `K40` [PPM].

Spatial discretization follows a halfwidth step of `dr` [μm], and temporal
discretization follows the age steps specified by the `agesteps` range,
in Ma.
"""
struct PlanarAr{T<:AbstractFloat} <: ArgonSample{T}
    age::T                      # [Ma] helium age
    age_sigma::T                # [Ma] helium age uncertainty (one-sigma)
    offset::T                   # [C] temperature offset relative to the surface
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
    isnothing(tsteps) && isnothing(agesteps) && @error "At least one of `tsteps` or `agesteps` is required"
    isnothing(tsteps) && (tsteps=reverse(agesteps))
    isnothing(agesteps) && (agesteps=reverse(tsteps))
    agesteps, tsteps = floatrange(agesteps), floatrange(tsteps)
    @assert issorted(tsteps)

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
        agesteps::AbstractRange,
        tsteps::AbstractRange=reverse(agesteps),
    )
    ```
    Construct a `MultipleDomain` diffusion chronometer given an observed argon
    release spectrum, degassing schedule, where domain is represented by a 
    `PlanarAr` or `SphericalAr` chronometer.

    Domain diffusivity and volume parameters must be supplied as vectors
    `Ea` [kJ/mol], `lnD0a2` [log(1/s)], and `volume_fraction` [unitless]
    obtained by separately fitting the release spectrum (the former two
    as an `MDDiffusivity` object).

    See also: `MDDiffusivity`, `PlanarAr`, `SphericalAr`, `degas!`
    """
    struct MultipleDomain{T<:AbstractFloat, C<:Union{SphericalAr{T}, PlanarAr{T}}} <: AbsoluteChronometer{T}
        age::Vector{T}                      # [Ma] measured ages at each degassing step
        age_sigma::Vector{T}                # [Ma] measured age uncertainties at each degassing step
        fraction_experimental::Vector{T}    # [unitless] fraction of total Ar-39 released each degassing step
        fraction_experimental_sigma::T      # [unitless] uncertainty in degassing fraction
        midpoint_experimental::Vector{T}    # [unitless] midpoint of fraction_experimental for each step
        tsteps_experimental::Vector{T}      # [s] time steps of experimental heating schedule
        Tsteps_experimental::Vector{T}      # [C] temperature steps of experimental heating schedule
        fit::BitVector                      # [Bool] Whether or not each step should be used in inversion
        offset::T                           # [C] temperature offset relative to the surface
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
            agesteps::AbstractRange,
            tsteps::AbstractRange=reverse(agesteps),
        )
        # Check input arrays are the right size and ordered properly
        @assert eachindex(age) == eachindex(age_sigma) == eachindex(fraction_experimental) == eachindex(tsteps_experimental) == eachindex(Tsteps_experimental)
        @assert issorted(tsteps_experimental, lt=<=) "Degassing time steps must be in strictly increasing order"
        @assert all(x->0<=x<=1, volume_fraction) "All \"fraction degassed\" values must be between 0 and 1"

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


## --- Internal functions to get values and uncertainties from any chronometers

val(x::AbsoluteChronometer{T}) where {T} = x.age::T
err(x::AbsoluteChronometer{T}) where {T} = x.age_sigma::T
val(x::MultipleDomain{T}) where {T} = nanmean(x.age, @.(x.fit/x.age_sigma^2))::T
err(x::MultipleDomain{T}) where {T} = nanstd(x.age, @.(x.fit/x.age_sigma^2))::T
val(x::FissionTrackLength{T}) where {T} = x.length::T
err(x::FissionTrackLength{T}) where {T} = zero(T)
val(x::ApatiteTrackLength{T}) where {T} = x.lcmod::T


## -- Functions related to age and age uncertinty of absolute chronometers

# Get age and age sigma from a vector of chronometers
function get_age(x::AbstractArray{<:Chronometer{T}}, ::Type{C}=AbsoluteChronometer{T}) where {T<:AbstractFloat, C<:AbsoluteChronometer}
    result = sizehint!(T[], length(x))
    for xᵢ in x
        if isa(xᵢ, C)
            push!(result, val(xᵢ))
        end
    end
    return result
end
function get_age_sigma(x::AbstractArray{<:Chronometer{T}}, ::Type{C}=AbsoluteChronometer{T}) where {T<:AbstractFloat, C<:AbsoluteChronometer}
    result = sizehint!(T[], length(x))
    for xᵢ in x
        if isa(xᵢ, C)
            push!(result, err(xᵢ))
        end
    end
    return result
end

function empiricaluncertainty!(σcalc::AbstractVector{T}, x::AbstractArray{<:Chronometer{T}}, ::Type{C}; fraction::Number=1/sqrt(2), sigma_eU::Number = ((C<:ZirconHe) ? 100. : 10.)) where {T<:AbstractFloat, C<:HeliumSample}
    @assert eachindex(σcalc) == eachindex(x)
    @assert 0 <= fraction <= 1
    t = isa.(x, C)
    eU_C = eU.(x[t])
    ages_C = get_age(x, C)
    for i ∈ findall(t)
        nearest_eU = minimum(j->abs(eU(x[j]) - eU(x[i])), setdiff(findall(t), i))
        W = normpdf.(eU(x[i]), max(sigma_eU, nearest_eU/2), eU_C)
        # Assume some fraction of weighted variance is from unknown external uncertainty
        σₑ = nanstd(ages_C, W) * fraction   # External uncertainty (est)
        σcalc[i] = max(σₑ, σcalc[i])
    end
    return x
end

## --- Importing of chronometers
"""
```julia
chronometers([T=Float64], data, model)
```
Construct a vector of `Chronometer` objects given a dataset `data`
and model parameters `model`
"""
chronometers(ds, model; kwargs...) = chronometers(Float64, ds, model; kwargs...)
function chronometers(T::Type{<:AbstractFloat}, ds, model;
        zirconvolumeweighting = :cylindrical,
        apatitevolumeweighting = :cylindrical,
    )
    agesteps = floatrange(model.agesteps)
    tsteps = floatrange(model.tsteps)
    dr = haskey(model, :dr) ? model.dr : one(T)
    @assert issorted(tsteps)
    @assert tsteps == reverse(agesteps)

    haskey(ds, :mineral) || @error "dataset must contain a column labeled `mineral`"
    mineral = ds.mineral
    crystage = if haskey(ds, :crystallization_age_Ma)        
        ds.crystallization_age_Ma
    elseif haskey(ds, :crystAge) # Legacy option
        ds.crystAge
    else
        @error "dataset must contain a column labeled `crystallization age [Ma]`"
    end
    @assert eachindex(mineral) == eachindex(crystage)

    # Default damage models for each mineral
    zdm = (haskey(model, :zdm) ? model.zdm : ZRDAAM())::ZirconHeliumModel{T}
    adm = (haskey(model, :adm) ? model.adm : RDAAM())::ApatiteHeliumModel{T}
    zftm = (haskey(model, :zftm) ? model.zftm : Yamada2007PC())::ZirconAnnealingModel{T}
    mftm = (haskey(model, :mftm) ? model.mftm : Jones2021FA())::MonaziteAnnealingModel{T}
    aftm = (haskey(model, :aftm) ? model.aftm : Ketcham2007FC())::ApatiteAnnealingModel{T}

    chrons = Chronometer[]
    damodels = Model[]
    for i in eachindex(mineral)
        first_index = 1 + round(Int,(maximum(agesteps) - crystage[i])/step(tsteps))
        mineral = lowercase(string(ds.mineral[i]))

        if mineral == "zircon"
            # Zircon helium
            if haskey(ds, :raw_He_age_Ma) && haskey(ds, :raw_He_age_sigma_Ma) && (0 < ds.raw_He_age_sigma_Ma[i]/ds.raw_He_age_Ma[i])
                # Modern format
                c = ZirconHe(T;
                    age = ds.raw_He_age_Ma[i], 
                    age_sigma = ds.raw_He_age_sigma_Ma[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    r = ds.halfwidth_um[i], 
                    dr = dr, 
                    U238 = (haskey(ds, :U238_ppm) && !isnan(ds.U238_ppm[i])) ? ds.U238_ppm[i] : 0,
                    Th232 = (haskey(ds, :Th232_ppm) && !isnan(ds.Th232_ppm[i])) ? ds.Th232_ppm[i] : 0,
                    Sm147 = (haskey(ds, :Sm147_ppm) && !isnan(ds.Sm147_ppm[i])) ? ds.Sm147_ppm[i] : 0,
                    U238_matrix = (haskey(ds, :U238_matrix_ppm) && !isnan(ds.U238_matrix_ppm[i])) ? ds.U238_matrix_ppm[i] : 0,
                    Th232_matrix = (haskey(ds, :Th232_matrix_ppm) && !isnan(ds.Th232_matrix_ppm[i])) ? ds.Th232_matrix_ppm[i] : 0,
                    Sm147_matrix = (haskey(ds, :Sm147_matrix_ppm) && !isnan(ds.Sm147_matrix_ppm[i])) ? ds.Sm147_matrix_ppm[i] : 0,
                    agesteps = agesteps[first_index:end],
                    volumeweighting = zirconvolumeweighting,
                )
                push!(chrons, c)
                push!(damodels, zdm)
            end
            # Zircon fission track
            if haskey(ds, :FT_age_Ma) && haskey(ds, :FT_age_sigma_Ma) &&  (0 < ds.FT_age_sigma_Ma[i]/ds.FT_age_Ma[i])
                c = ZirconFT(T;
                    age = ds.FT_age_Ma[i], 
                    age_sigma = ds.FT_age_sigma_Ma[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    agesteps = agesteps[first_index:end],
                )
                push!(chrons, c)
                push!(damodels, zftm)
            end
            # Zircon fission track length
            if haskey(ds, :track_length_um) && (0 < ds.track_length_um[i])
                c = ZirconTrackLength(T;
                    length = ds.track_length_um[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    agesteps = agesteps[first_index:end],
                )
                push!(chrons, c)
                push!(damodels, zftm)
            end

        elseif mineral == "monazite"
            # Monazite fission track
            if haskey(ds, :FT_age_Ma) && haskey(ds, :FT_age_sigma_Ma) && (0 < ds.FT_age_sigma_Ma[i]/ds.FT_age_Ma[i])
                c = MonaziteFT(T;
                    age = ds.FT_age_Ma[i], 
                    age_sigma = ds.FT_age_sigma_Ma[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    agesteps = agesteps[first_index:end],
                )
                push!(chrons, c)
                push!(damodels, mftm)
            end
            # Monazite fission track length
            if haskey(ds, :track_length_um) && (0 < ds.track_length_um[i])
                c = MonaziteTrackLength(T;
                    length = ds.track_length_um[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    agesteps = agesteps[first_index:end],
                )
                push!(chrons, c)
                push!(damodels, mftm)
            end
            
        elseif mineral == "apatite"
            # Apatite helium
            if haskey(ds, :raw_He_age_Ma) && haskey(ds, :raw_He_age_sigma_Ma) && (0 < ds.raw_He_age_sigma_Ma[i]/ds.raw_He_age_Ma[i])
                # Modern format
                c = ApatiteHe(T;
                    age = ds.raw_He_age_Ma[i], 
                    age_sigma = ds.raw_He_age_sigma_Ma[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    r = ds.halfwidth_um[i], 
                    dr = dr, 
                    U238 = (haskey(ds, :U238_ppm) && !isnan(ds.U238_ppm[i])) ? ds.U238_ppm[i] : 0,
                    Th232 = (haskey(ds, :Th232_ppm) && !isnan(ds.Th232_ppm[i])) ? ds.Th232_ppm[i] : 0,
                    Sm147 = (haskey(ds, :Sm147_ppm) && !isnan(ds.Sm147_ppm[i])) ? ds.Sm147_ppm[i] : 0,
                    U238_matrix = (haskey(ds, :U238_matrix_ppm) && !isnan(ds.U238_matrix_ppm[i])) ? ds.U238_matrix_ppm[i] : 0,
                    Th232_matrix = (haskey(ds, :Th232_matrix_ppm) && !isnan(ds.Th232_matrix_ppm[i])) ? ds.Th232_matrix_ppm[i] : 0,
                    Sm147_matrix = (haskey(ds, :Sm147_matrix_ppm) && !isnan(ds.Sm147_matrix_ppm[i])) ? ds.Sm147_matrix_ppm[i] : 0,
                    agesteps = agesteps[first_index:end],
                    volumeweighting = apatitevolumeweighting,
                )
                push!(chrons, c)
                push!(damodels, adm)
            end
            # Apatite fission track
            if haskey(ds, :FT_age_Ma) && haskey(ds, :FT_age_sigma_Ma) && (0 < ds.FT_age_sigma_Ma[i]/ds.FT_age_Ma[i])
                c = ApatiteFT(T;
                    age = ds.FT_age_Ma[i], 
                    age_sigma = ds.FT_age_sigma_Ma[i], 
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    dpar = haskey(ds, :dpar_um) ? ds.dpar_um[i] : NaN,
                    F = haskey(ds, :F_apfu) ? ds.F_apfu[i] : NaN,
                    Cl = haskey(ds, :Cl_apfu) ? ds.Cl_apfu[i] : NaN,
                    OH = haskey(ds, :OH_apfu) ? ds.OH_apfu[i] : NaN,
                    rmr0 =haskey(ds, :rmr0) ? ds.rmr0[i] : NaN,
                    agesteps = agesteps[first_index:end],
                )
                push!(chrons, c)
                push!(damodels, aftm)
            end
            # Apatite fission track length
            if haskey(ds, :track_length_um) && (0 < ds.track_length_um[i])
                c = ApatiteTrackLength(T;
                    length = ds.track_length_um[i], 
                    angle = (haskey(ds, :track_angle_degrees) && !isnan(ds.track_angle_degrees[i])) ? ds.track_angle_degrees[i] : 0,
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    dpar = haskey(ds, :dpar_um) ? ds.dpar_um[i] : NaN,
                    F = haskey(ds, :F_apfu) ? ds.F_apfu[i] : NaN,
                    Cl = haskey(ds, :Cl_apfu) ? ds.Cl_apfu[i] : NaN,
                    OH = haskey(ds, :OH_apfu) ? ds.OH_apfu[i] : NaN,
                    rmr0 =haskey(ds, :rmr0) ? ds.rmr0[i] : NaN,
                    agesteps = agesteps[first_index:end],
                )
                push!(chrons, c)
                push!(damodels, aftm)
            end
        elseif haskey(ds, :D0_cm_2_s) && haskey(ds, :Ea_kJ_mol) && (0 < ds.D0_cm_2_s[i]) && (0 < ds.Ea_kJ_mol[i])
            geometry = haskey(ds, :geometry) ? lowercase(string(ds.geometry[i])) : "spherical"
            if (geometry == "slab") || (geometry == "planar")
                # Planar slab helium
                if haskey(ds, :raw_He_age_Ma) && haskey(ds, :raw_He_age_sigma_Ma) && (0 < ds.raw_He_age_sigma_Ma[i]/ds.raw_He_age_Ma[i])
                    c = PlanarHe(T;
                        age = ds.raw_He_age_Ma[i], 
                        age_sigma = ds.raw_He_age_sigma_Ma[i], 
                        offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                        stoppingpower = alphastoppingpower(ds.mineral[i]),
                        r = ds.halfwidth_um[i], 
                        dr = dr, 
                        U238 = (haskey(ds, :U238_ppm) && !isnan(ds.U238_ppm[i])) ? ds.U238_ppm[i] : 0,
                        Th232 = (haskey(ds, :Th232_ppm) && !isnan(ds.Th232_ppm[i])) ? ds.Th232_ppm[i] : 0,
                        Sm147 = (haskey(ds, :Sm147_ppm) && !isnan(ds.Sm147_ppm[i])) ? ds.Sm147_ppm[i] : 0,
                        U238_matrix = (haskey(ds, :U238_matrix_ppm) && !isnan(ds.U238_matrix_ppm[i])) ? ds.U238_matrix_ppm[i] : 0,
                        Th232_matrix = (haskey(ds, :Th232_matrix_ppm) && !isnan(ds.Th232_matrix_ppm[i])) ? ds.Th232_matrix_ppm[i] : 0,
                        Sm147_matrix = (haskey(ds, :Sm147_matrix_ppm) && !isnan(ds.Sm147_matrix_ppm[i])) ? ds.Sm147_matrix_ppm[i] : 0,
                        agesteps = agesteps[first_index:end],
                    )
                    dm = Diffusivity(
                        D0 = T(ds.D0_cm_2_s[i]),
                        D0_logsigma = T((haskey(ds, :D0_logsigma) && !isnan(ds.D0_logsigma[i])) ? ds.D0_logsigma[i] : log(2)/2),
                        Ea = T(ds.Ea_kJ_mol[i]),
                        Ea_logsigma = T((haskey(ds, :Ea_logsigma) && !isnan(ds.Ea_logsigma[i])) ? ds.Ea_logsigma[i] : log(2)/4),
                    )
                    push!(chrons, c)
                    push!(damodels, dm)
                end
                # Planar slab argon
                if haskey(ds, :raw_Ar_age_Ma) && haskey(ds, :raw_Ar_age_sigma_Ma) && (0 < ds.raw_Ar_age_sigma_Ma[i]/ds.raw_Ar_age_Ma[i])
                    c = PlanarAr(T;
                        age = ds.raw_Ar_age_Ma[i], 
                        age_sigma = ds.raw_Ar_age_sigma_Ma[i], 
                        offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                        r = ds.halfwidth_um[i], 
                        dr = dr, 
                        K40 = (haskey(ds, :K40_ppm) && !isnan(ds.K40_ppm[i])) ? ds.K40_ppm[i] : 16.34,
                        agesteps = agesteps[first_index:end],
                    )
                    dm = Diffusivity(
                        D0 = T(ds.D0_cm_2_s[i]),
                        D0_logsigma = T((haskey(ds, :D0_logsigma) && !isnan(ds.D0_logsigma[i])) ? ds.D0_logsigma[i] : log(2)/2),
                        Ea = T(ds.Ea_kJ_mol[i]),
                        Ea_logsigma = T((haskey(ds, :Ea_logsigma) && !isnan(ds.Ea_logsigma[i])) ? ds.Ea_logsigma[i] : log(2)/4),
                    )
                    push!(chrons, c)
                    push!(damodels, dm)
                end
            else
                (geometry === "sphere") || (geometry === "spherical") || @warn "Geometry \"$geometry\" not recognized in row $i, defaulting to spherical"
                # Spherical helium
                if haskey(ds, :raw_He_age_Ma) && haskey(ds, :raw_He_age_sigma_Ma) && (0 < ds.raw_He_age_sigma_Ma[i]/ds.raw_He_age_Ma[i])
                    c = SphericalHe(T;
                        age = ds.raw_He_age_Ma[i], 
                        age_sigma = ds.raw_He_age_sigma_Ma[i], 
                        offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                        stoppingpower = alphastoppingpower(ds.mineral[i]),
                        r = ds.halfwidth_um[i], 
                        dr = dr, 
                        U238 = (haskey(ds, :U238_ppm) && !isnan(ds.U238_ppm[i])) ? ds.U238_ppm[i] : 0,
                        Th232 = (haskey(ds, :Th232_ppm) && !isnan(ds.Th232_ppm[i])) ? ds.Th232_ppm[i] : 0,
                        Sm147 = (haskey(ds, :Sm147_ppm) && !isnan(ds.Sm147_ppm[i])) ? ds.Sm147_ppm[i] : 0,
                        U238_matrix = (haskey(ds, :U238_matrix_ppm) && !isnan(ds.U238_matrix_ppm[i])) ? ds.U238_matrix_ppm[i] : 0,
                        Th232_matrix = (haskey(ds, :Th232_matrix_ppm) && !isnan(ds.Th232_matrix_ppm[i])) ? ds.Th232_matrix_ppm[i] : 0,
                        Sm147_matrix = (haskey(ds, :Sm147_matrix_ppm) && !isnan(ds.Sm147_matrix_ppm[i])) ? ds.Sm147_matrix_ppm[i] : 0,
                        agesteps = agesteps[first_index:end],
                    )
                    dm = Diffusivity(
                        D0 = T(ds.D0_cm_2_s[i]),
                        D0_logsigma = T((haskey(ds, :D0_logsigma) && !isnan(ds.D0_logsigma[i])) ? ds.D0_logsigma[i] : log(2)/2),
                        Ea = T(ds.Ea_kJ_mol[i]),
                        Ea_logsigma = T((haskey(ds, :Ea_logsigma) && !isnan(ds.Ea_logsigma[i])) ? ds.Ea_logsigma[i] : log(2)/4),
                    )
                    push!(chrons, c)
                    push!(damodels, dm)
                end
                # Spherical argon
                if haskey(ds, :raw_Ar_age_Ma) && haskey(ds, :raw_Ar_age_sigma_Ma) && (0 < ds.raw_Ar_age_sigma_Ma[i]/ds.raw_Ar_age_Ma[i])
                    c = SphericalAr(T;
                        age = ds.raw_Ar_age_Ma[i], 
                        age_sigma = ds.raw_Ar_age_sigma_Ma[i], 
                        offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                        r = ds.halfwidth_um[i], 
                        dr = dr, 
                        K40 = (haskey(ds, :K40_ppm) && !isnan(ds.K40_ppm[i])) ? ds.K40_ppm[i] : 16.34,
                        agesteps = agesteps[first_index:end],
                    )
                    dm = Diffusivity(
                        D0 = T(ds.D0_cm_2_s[i]),
                        D0_logsigma = T((haskey(ds, :D0_logsigma) && !isnan(ds.D0_logsigma[i])) ? ds.D0_logsigma[i] : log(2)/2),
                        Ea = T(ds.Ea_kJ_mol[i]),
                        Ea_logsigma = T((haskey(ds, :Ea_logsigma) && !isnan(ds.Ea_logsigma[i])) ? ds.Ea_logsigma[i] : log(2)/4),
                    )
                    push!(chrons, c)
                    push!(damodels, dm)
                end
            end
        elseif haskey(ds, :mdd_file) && !isempty(ds.mdd_file[i])
            mdds = importdataset(ds.mdd_file[i], importas=:Tuple)
            geometry = haskey(ds, :geometry) ? lowercase(string(ds.geometry[i])) : "spherical"
            r = (haskey(ds, :halfwidth_um) && !isnan(ds.halfwidth_um[i])) ? ds.halfwidth_um[i] : 100
            if (geometry == "slab") || (geometry == "planar")
                c = MultipleDomain(T, PlanarAr;
                    age = mdds.age_Ma,
                    age_sigma = mdds.age_sigma_Ma,
                    fraction_experimental = mdds.fraction_degassed,
                    tsteps_experimental = issorted(mdds.time_s, lt=<=) ? mdds.time_s : cumsum(mdds.time_s),
                    Tsteps_experimental = mdds.temperature_C,
                    fit = mdds.fit,
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    r = r,
                    dr = dr, 
                    volume_fraction = mdds.volume_fraction[.!isnan.(mdds.volume_fraction)],
                    K40 = (haskey(ds, :K40_ppm) && !isnan(ds.K40_ppm[i])) ? ds.K40_ppm[i] : 16.34,
                    agesteps = agesteps[first_index:end],
                )
            else
                (geometry === "spherical") || @warn "Geometry \"$geometry\" not recognized in row $i, defaulting to spherical"
                c = MultipleDomain(T, SphericalAr;
                    age = mdds.age_Ma,
                    age_sigma = mdds.age_sigma_Ma,
                    fraction_experimental = mdds.fraction_degassed,
                    tsteps_experimental = issorted(mdds.time_s, lt=<=) ? mdds.time_s : cumsum(mdds.time_s),
                    Tsteps_experimental = mdds.temperature_C,
                    fit = mdds.fit,
                    offset = (haskey(ds, :offset_C) && !isnan(ds.offset_C[i])) ? ds.offset_C[i] : 0,
                    r = r,
                    dr = dr, 
                    volume_fraction = mdds.volume_fraction[.!isnan.(mdds.volume_fraction)],
                    K40 = (haskey(ds, :K40_ppm) && !isnan(ds.K40_ppm[i])) ? ds.K40_ppm[i] : 16.34,
                    agesteps = agesteps[first_index:end],
                )
            end
            tdomains = .!isnan.(mdds.lnD0_a_2)
            dm = MDDiffusivity(
                D0 = (T.(exp.(mdds.lnD0_a_2[tdomains]).*(r/10000)^2)...,),
                D0_logsigma = (T.(haskey(mdds, :lnD0_a_2_sigma) ? mdds.lnD0_a_2_sigma[tdomains] : fill(log(2)/2, count(tdomains)))...,),
                Ea = (T.(mdds.Ea_kJ_mol[tdomains])...,),
                Ea_logsigma = (T.(haskey(mdds, :Ea_logsigma) ? mdds.Ea_logsigma[tdomains] : fill(log(2)/4, count(tdomains)))...,),
            )
            push!(chrons, c)
            push!(damodels, dm)
        end
    end

    isempty(chrons) && @error "No chronometers found"
    return unionize(chrons), unionize(damodels)
end
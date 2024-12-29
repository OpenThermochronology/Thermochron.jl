# Abstract type to include any number of mineral chronometers (zircon, apatite, etc.)
abstract type Chronometer{T} end

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
eU(x::Chronometer{T}) where {T<:AbstractFloat} = T(NaN)
function eU(x::HeliumSample{T}) where {T<:AbstractFloat}
    # Convert from atoms/g to ppm
    eu = nanmean(x.r238U) / (6.022E23 / 1E6 / 238)
    eu += 0.238*nanmean(x.r232Th) / (6.022E23 / 1E6 / 232)
    eu += 0.012*nanmean(x.r147Sm) / (6.022E23 / 1E6 / 147)
    return T(eu)
end

## --- Internal functions to get values and uncertainties from any chronometers
val(x::AbsoluteChronometer{T}) where {T} = x.age::T
err(x::AbsoluteChronometer{T}) where {T} = x.age_sigma::T
val(x::FissionTrackLength{T}) where {T} = x.lcmod::T
err(x::FissionTrackLength{T}) where {T} = zero(T)


## --- Fission track sample types

struct ApatiteTrackLength{T<:AbstractFloat} <: FissionTrackLength{T}
    length::T               # [um]
    angle::T                # [degrees]
    lcmod::T                # [um]
    agesteps::FloatRange    # [Ma]
    tsteps::FloatRange      # [Ma]
    r::Vector{T}            # [unitless]
    pr::Vector{T}           # [unitless]
    dpar::T                 # [μm]
    F::T                    # [APFU]
    Cl::T                   # [APFU]
    OH::T                   # [APFU]
    rmr0::T                 # [unitless]
end
function ApatiteTrackLength(T::Type{<:AbstractFloat}=Float64; 
    length=T(NaN), 
    angle=T(NaN), 
    lcmod=lcmod(length, angle),
    agesteps, 
    tsteps=reverse(agesteps), 
    r=zeros(T, size(agesteps)),
    pr=zeros(T, size(agesteps)),
    dpar=T(NaN), 
    F=T(NaN), 
    Cl=T(NaN), 
    OH=T(NaN), 
    rmr0=T(NaN),
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
ApatiteTrackLength(
    T(length),
    T(angle),
    T(lcmod),
    floatrange(agesteps),
    floatrange(tsteps),
    T.(r),
    T.(pr),
    T(dpar),
    T(F),
    T(Cl),
    T(OH),
    T(rmr0),
)
end

struct ApatiteFT{T<:AbstractFloat} <: FissionTrackSample{T}
    age::T                  # [Ma]
    age_sigma::T            # [Ma]
    agesteps::FloatRange    # [Ma]
    tsteps::FloatRange      # [Ma]
    dpar::T                 # [μm]
    F::T                    # [APFU]
    Cl::T                   # [APFU]
    OH::T                   # [APFU]
    rmr0::T                 # [unitless]
end
function ApatiteFT(T::Type{<:AbstractFloat}=Float64; 
        age=T(NaN), 
        age_sigma=T(NaN), 
        agesteps, 
        tsteps=reverse(agesteps), 
        dpar=T(NaN),
        F=T(NaN),
        Cl=T(NaN),
        OH=T(NaN),
        rmr0=T(NaN),
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
        floatrange(agesteps),
        floatrange(tsteps),
        T(dpar),
        T(F),
        T(Cl),
        T(OH),
        T(rmr0),
    )
end

## --- Helium sample types
"""
```julia
ZirconHe(T=Float64;
    age::Number=T(NaN),
    age_sigma::Number=T(NaN),
    r::Number=one(T),
    dr::Number, 
    U238::Number,
    Th232::Number,
    Sm147::Number=zero(T),
    agesteps::AbstractRange
)
```
Construct a `ZirconHe` chronometer representing a zircon with a raw 
helium age of `age` ± `age_sigma` Ma, a  radius of `r` μm, and uniform 
U, Th and Sm concentrations in ppm specified by `U238`, `Th232`, and `Sm147`. 
A present day U-235/U-238 ratio of 1/137.818 is assumed.

Spatial discretization follows a radius step of `dr` μm, and temporal
discretization follows the age steps specified by the `agesteps` range,
in Ma.
"""
# Concretely-typed immutable struct to hold information about a single zircon (Helium) crystal
struct ZirconHe{T<:AbstractFloat} <: HeliumSample{T}
    age::T
    age_sigma::T
    agesteps::FloatRange
    tsteps::FloatRange
    rsteps::FloatRange
    redges::FloatRange
    nrsteps::Int
    r238U::Vector{T}
    r235U::Vector{T}
    r232Th::Vector{T}
    r147Sm::Vector{T}
    alphadeposition::Matrix{T}
    alphadamage::Matrix{T}
    pr::Matrix{T}
    annealeddamage::Matrix{T}
    u::Matrix{T}
    β::Vector{T}
    Dz::Vector{T}
    DN17::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{T, Tridiagonal{T, Vector{T}}, Vector{Int64}}
    y::Vector{T}
end
# Constructor for the ZirconHe type, given grain radius, U and Th concentrations and t-T discretization information
function ZirconHe(T::Type{<:AbstractFloat}=Float64;
        age::Number=T(NaN),
        age_sigma::Number=T(NaN),
        r::Number=one(T),
        dr::Number, 
        U238::Number,
        Th232::Number,
        Sm147::Number=zero(T),
        agesteps::AbstractRange
    )
    
    # Temporal discretization
    agesteps = floatrange(agesteps)
    tsteps = reverse(agesteps)
    @assert issorted(tsteps)

    # crystal size and spatial discretization
    rsteps = floatrange(0+dr/2 : dr : r-dr/2)
    redges = floatrange(     0 : dr : r     )   # Edges of each radius element
    nrsteps = length(rsteps)+2                  # number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    relvolumes = (redges[2:end].^3 - redges[1:end-1].^3)/redges[end]^3 # Relative volume fraction of spherical shell corresponding to each radius element

    # Alpha stopping distances for each isotope in each decay chain, from
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alpharadii238U = (11.78, 14.09, 13.73, 14.13, 17.32, 16.69, 28.56, 16.48,)
    alpharadii235U = (12.58, 15.04, 19.36, 18.06, 23.07, 26.87, 22.47,)
    alpharadii232Th = (10.99, 16.67, 18.16, 17.32, 23.61, 29.19,)
    # Ketchem et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alpharadii147Sm = (4.76,)


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

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance
    dint = zeros(T, length(redges) - 1)
    
    r238UHe = zeros(T, size(r238U))
    @inbounds for ri in eachindex(rsteps, relvolumes, r238U)
        for i in eachindex(alpharadii238U)
            # Effective radial alpha deposition from 238U
            intersectiondensity!(dint,redges,relvolumes,alpharadii238U[i],rsteps[ri])
            @. r238UHe += relvolumes[ri] * dint * r238U[ri]
        end
    end

    # Effective radial alpha deposition from U-235
    r235UHe = zeros(T, size(r235U))
    @inbounds for ri in eachindex(rsteps, relvolumes, r235U)
        for i in eachindex(alpharadii235U)
            # Effective radial alpha deposition from 235U
            intersectiondensity!(dint, redges,relvolumes,alpharadii235U[i],rsteps[ri])
            @. r235UHe += relvolumes[ri] * dint * r235U[ri]
        end
    end

    # Effective radial alpha deposition from Th-232
    r232ThHe = zeros(T, size(r232Th))
    @inbounds for ri in eachindex(rsteps, relvolumes, r232Th)
        for i in eachindex(alpharadii232Th)
            # Effective radial alpha deposition from 232Th
            intersectiondensity!(dint, redges,relvolumes,alpharadii232Th[i],rsteps[ri])
            @. r232ThHe += relvolumes[ri] * dint * r232Th[ri]
        end
    end

    # Effective radial alpha deposition from Sm-147
    r147SmHe = zeros(T, size(r147Sm))
    @inbounds for ri in eachindex(rsteps, relvolumes, r147Sm)
        for i in eachindex(alpharadii147Sm)
            intersectiondensity!(dint, redges,relvolumes,alpharadii147Sm[i],rsteps[ri])
            @. r147SmHe += relvolumes[ri] * dint * r147Sm[ri]
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

    # Allocate additional variables that will be needed for Crank-Nicholson
    annealeddamage = similar(alphadamage)
    β = zeros(T, nrsteps) # First row of annealeddamage

    # Allocate arrays for diffusivities
    Dz = zeros(T, length(tsteps))
    DN17 = zeros(T, length(tsteps))

    # Allocate output matrix for all timesteps
    u = zeros(T, nrsteps, length(tsteps))

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

    return ZirconHe(
        T(age),
        T(age_sigma),
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
    age::Number=T(NaN),
    age_sigma::Number=T(NaN),
    r::Number, 
    dr::Number=one(T), 
    U238::Number, 
    Th232::Number, 
    Sm147::Number=zero(T), 
    agesteps::AbstractRange,
)
```
Construct an `ApatiteHe` chronometer representing an apatite with a raw 
helium age of `age` ± `age_sigma` Ma, a  radius of `r` μm, and uniform 
U, Th and Sm concentrations in ppm specified by `U238`, `Th232`, and `Sm147`. 
A present day U-235/U-238 ratio of 1/137.818 is assumed.

Spatial discretization follows a radius step of `dr` μm, and temporal
discretization follows the age steps specified by the `agesteps` range,
in Ma.
"""
# Concretely-typed immutable struct to hold information about a single apatite (Helium) crystal
struct ApatiteHe{T<:AbstractFloat} <: HeliumSample{T}
    age::T
    age_sigma::T
    agesteps::FloatRange
    tsteps::FloatRange
    rsteps::FloatRange
    redges::FloatRange
    nrsteps::Int
    r238U::Vector{T}
    r235U::Vector{T}
    r232Th::Vector{T}
    r147Sm::Vector{T}
    alphadeposition::Matrix{T}
    alphadamage::Matrix{T}
    pr::Matrix{T}
    annealeddamage::Matrix{T}
    u::Matrix{T}
    β::Vector{T}
    DL::Vector{T}
    Dtrap::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{T, Tridiagonal{T, Vector{T}}, Vector{Int64}}
    y::Vector{T}
end
# Constructor for the ApatiteHe type, given grain radius, U and Th concentrations and t-T discretization information
function ApatiteHe(T::Type{<:AbstractFloat}=Float64;
        age::Number=T(NaN),
        age_sigma::Number=T(NaN),
        r::Number, 
        dr::Number=one(T), 
        U238::Number, 
        Th232::Number, 
        Sm147::Number=zero(T), 
        agesteps::AbstractRange,
    )

    # Temporal discretization
    agesteps = floatrange(agesteps)
    tsteps = reverse(agesteps)
    @assert issorted(tsteps)

    # crystal size and spatial discretization
    rsteps = floatrange(0+dr/2 : dr : r-dr/2)
    redges = floatrange(     0 : dr : r     )   # Edges of each radius element
    nrsteps = length(rsteps)+2                  # number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    relvolumes = (redges[2:end].^3 - redges[1:end-1].^3)/redges[end]^3 # Relative volume fraction of spherical shell corresponding to each radius element

    # Alpha stopping distances for each isotope in each decay chain, from
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alpharadii238U = (13.54, 16.26, 15.84, 16.31, 20.09, 22.89, 33.39, 19.10,)
    alpharadii235U = (14.48, 17.39, 22.5, 20.97, 26.89, 31.40, 26.18,)
    alpharadii232Th = (12.60, 19.32, 21.08, 20.09, 27.53, 34.14,)
    # Ketchem et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alpharadii147Sm = (5.93,)

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

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance
    dint = zeros(T, length(redges) - 1)

    # Effective radial alpha deposition from U-238
    r238UHe = zeros(T, size(r238U))
    @inbounds for ri in eachindex(rsteps, relvolumes, r238U)
        for i in eachindex(alpharadii238U)
            # Effective radial alpha deposition from 238U
            intersectiondensity!(dint,redges,relvolumes,alpharadii238U[i],rsteps[ri])
            @. r238UHe += relvolumes[ri] * dint * r238U[ri]
        end
    end

    # Effective radial alpha deposition from U-235
    r235UHe = zeros(T, size(r235U))
    @inbounds for ri in eachindex(rsteps, relvolumes, r235U)
        for i in eachindex(alpharadii235U)
            # Effective radial alpha deposition from 235U
            intersectiondensity!(dint, redges,relvolumes,alpharadii235U[i],rsteps[ri])
            @. r235UHe += relvolumes[ri] * dint * r235U[ri]
        end
    end

    # Effective radial alpha deposition from Th-232
    r232ThHe = zeros(T, size(r232Th))
    @inbounds for ri in eachindex(rsteps, relvolumes, r232Th)
        for i in eachindex(alpharadii232Th)
            # Effective radial alpha deposition from 232Th
            intersectiondensity!(dint, redges,relvolumes,alpharadii232Th[i],rsteps[ri])
            @. r232ThHe += relvolumes[ri] * dint * r232Th[ri]
        end
    end

    # Effective radial alpha deposition from Sm-147
    r147SmHe = zeros(T, size(r147Sm))
    @inbounds for ri in eachindex(rsteps, relvolumes, r147Sm)
        for i in eachindex(alpharadii147Sm)
            intersectiondensity!(dint, redges,relvolumes,alpharadii147Sm[i],rsteps[ri])
            @. r147SmHe += relvolumes[ri] * dint * r147Sm[ri]
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

    # Allocate additional variables that will be needed for Crank-Nicholson
    annealeddamage = similar(alphadamage)
    β = zeros(T, nrsteps) # First row of annealeddamage

    # Allocate arrays for diffusivities
    DL = zeros(T, length(tsteps))
    Dtrap = zeros(T, length(tsteps))

    # Allocate output matrix for all timesteps
    u = zeros(T, nrsteps, length(tsteps))

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

    return ApatiteHe(
        T(age),
        T(age_sigma),
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
GenericHe(T=Float64;
    age::Number=T(NaN),
    age_sigma::Number=T(NaN),
    D0::Number,
    Ea::Number,
    r::Number, 
    dr::Number=one(T), 
    U238::Number, 
    Th232::Number, 
    Sm147::Number=zero(T), 
    agesteps::AbstractRange,
)
```
Construct an `GenericHe` chronometer representing a mineral with a raw 
helium age of `age` ± `age_sigma` Ma, a uniform diffusivity specified by
a D₀ of `D0` cm^2/sec and an activation energy of `Ea` kJ/mol, a 
radius of `r` μm, and uniform U, Th and Sm concentrations in ppm
specified by `U238`, `Th232`, and `Sm147`. (A present day U-235/U-238 
ratio of 1/137.818 is assumed)

Spatial discretization follows a radius step of `dr` μm, and temporal
discretization follows the age steps specified by the `agesteps` range,
in Ma.
"""
# Concretely-typed immutable struct to hold information about a single mineral (Helium) crystal
struct GenericHe{T<:AbstractFloat} <: HeliumSample{T}
    age::T
    age_sigma::T
    D0::T
    Ea::T
    agesteps::FloatRange
    tsteps::FloatRange
    rsteps::FloatRange
    redges::FloatRange
    nrsteps::Int
    r238U::Vector{T}
    r235U::Vector{T}
    r232Th::Vector{T}
    r147Sm::Vector{T}
    alphadeposition::Matrix{T}
    u::Matrix{T}
    β::Vector{T}
    De::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{T, Tridiagonal{T, Vector{T}}, Vector{Int64}}
    y::Vector{T}
end
# Constructor for the GenericHe type, given grain radius, U and Th concentrations and t-T discretization information
function GenericHe(T::Type{<:AbstractFloat}=Float64;
        age::Number=T(NaN),
        age_sigma::Number=T(NaN),
        D0::Number,
        Ea::Number,
        r::Number, 
        dr::Number=one(T), 
        U238::Number, 
        Th232::Number, 
        Sm147::Number=zero(T), 
        agesteps::AbstractRange,
    )

    # Temporal discretization
    agesteps = floatrange(agesteps)
    tsteps = reverse(agesteps)
    @assert issorted(tsteps)

    # crystal size and spatial discretization
    rsteps = floatrange(0+dr/2 : dr : r-dr/2)
    redges = floatrange(     0 : dr : r     )   # Edges of each radius element
    nrsteps = length(rsteps)+2                  # number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    relvolumes = (redges[2:end].^3 - redges[1:end-1].^3)/redges[end]^3 # Relative volume fraction of spherical shell corresponding to each radius element

    # Alpha stopping distances for each isotope in each decay chain, from
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alpharadii238U = (13.54, 16.26, 15.84, 16.31, 20.09, 22.89, 33.39, 19.10,)
    alpharadii235U = (14.48, 17.39, 22.5, 20.97, 26.89, 31.40, 26.18,)
    alpharadii232Th = (12.60, 19.32, 21.08, 20.09, 27.53, 34.14,)
    # Ketchem et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alpharadii147Sm = (5.93,)

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

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance
    dint = zeros(T, length(redges) - 1)

    # Effective radial alpha deposition from U-238
    r238UHe = zeros(T, size(r238U))
    @inbounds for ri in eachindex(rsteps, relvolumes, r238U)
        for i in eachindex(alpharadii238U)
            # Effective radial alpha deposition from 238U
            intersectiondensity!(dint,redges,relvolumes,alpharadii238U[i],rsteps[ri])
            @. r238UHe += relvolumes[ri] * dint * r238U[ri]
        end
    end

    # Effective radial alpha deposition from U-235
    r235UHe = zeros(T, size(r235U))
    @inbounds for ri in eachindex(rsteps, relvolumes, r235U)
        for i in eachindex(alpharadii235U)
            # Effective radial alpha deposition from 235U
            intersectiondensity!(dint, redges,relvolumes,alpharadii235U[i],rsteps[ri])
            @. r235UHe += relvolumes[ri] * dint * r235U[ri]
        end
    end

    # Effective radial alpha deposition from Th-232
    r232ThHe = zeros(T, size(r232Th))
    @inbounds for ri in eachindex(rsteps, relvolumes, r232Th)
        for i in eachindex(alpharadii232Th)
            # Effective radial alpha deposition from 232Th
            intersectiondensity!(dint, redges,relvolumes,alpharadii232Th[i],rsteps[ri])
            @. r232ThHe += relvolumes[ri] * dint * r232Th[ri]
        end
    end

    # Effective radial alpha deposition from Sm-147
    r147SmHe = zeros(T, size(r147Sm))
    @inbounds for ri in eachindex(rsteps, relvolumes, r147Sm)
        for i in eachindex(alpharadii147Sm)
            intersectiondensity!(dint, redges,relvolumes,alpharadii147Sm[i],rsteps[ri])
            @. r147SmHe += relvolumes[ri] * dint * r147Sm[ri]
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

    # Allocate additional variables that will be needed for Crank-Nicholson
    β = zeros(T, nrsteps) # First row of annealeddamage

    # Allocate arrays for diffusivities
    De = zeros(T, length(tsteps))

    # Allocate output matrix for all timesteps
    u = zeros(T, nrsteps, length(tsteps))

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

    return GenericHe(
        T(age),
        T(age_sigma),
        T(D0),
        T(Ea),
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

## -- Functions related to age and age uncertinty of absolute chronometers

# Get age and age sigma from a vector of chronometers
function get_age(x::AbstractArray{<:Chronometer{T}}, ::Type{C}=AbsoluteChronometer{T}) where {T<:AbstractFloat, C<:AbsoluteChronometer}
    result = sizehint!(T[], length(x))
    for xᵢ in x
        if isa(xᵢ, C)
            push!(result, xᵢ.age)
        end
    end
    return result
end
function get_age_sigma(x::AbstractArray{<:Chronometer{T}}, ::Type{C}=AbsoluteChronometer{T}) where {T<:AbstractFloat, C<:AbsoluteChronometer}
    result = sizehint!(T[], length(x))
    for xᵢ in x
        if isa(xᵢ, C)
            push!(result, xᵢ.age_sigma)
        end
    end
    return result
end
# Set age and age sigma from a vector of chronometers
function set_age!(x::AbstractArray{<:Chronometer{T}}, ages::AbstractArray{T}, ::Type{C}=AbsoluteChronometer{T}) where {T<:AbstractFloat, C<:AbsoluteChronometer}
    @assert count(xᵢ->isa(xᵢ, C), x) == length(ages)
    iₐ = firstindex(ages)
    for iₓ in eachindex(x)
        if isa(x[iₓ], C)
            x[iₓ] = setproperty!!(x[iₓ], :age, ages[iₐ])
            iₐ += 1
        end
    end
    return x
end
function set_age_sigma!(x::AbstractArray{<:Chronometer{T}}, age_sigmas::AbstractArray{T}, ::Type{C}=AbsoluteChronometer{T}) where {T<:AbstractFloat, C<:AbsoluteChronometer}
    @assert count(xᵢ->isa(xᵢ, C), x) == length(age_sigmas)
    iₐ = firstindex(age_sigmas)
    for iₓ in eachindex(x)
        if isa(x[iₓ], C)
            x[iₓ] = setproperty!!(x[iₓ], :age_sigma, age_sigmas[iₐ])
            iₐ += 1
        end
    end
    return x
end

function empiricaluncertainty!(x::AbstractArray{<:Chronometer{T}}, ::Type{C}; fraction::Number=1/sqrt(2), sigma_eU::Number = ((C<:ZirconHe) ? 100. : 10.)) where {T<:AbstractFloat, C<:HeliumSample}
    @assert 0 <= fraction <= 1
    t = isa.(x, C)
    eU_C = eU.(x[t])
    ages_C = get_age(x, C)
    for i ∈ findall(t)
        nearest_eU = minimum(j->abs(eU(x[j]) - eU(x[i])), setdiff(findall(t), i))
        W = normpdf.(eU(x[i]), max(sigma_eU, nearest_eU/2), eU_C)
        # Assume some fraction of weighted variance is from unknown external uncertainty
        σₑ = nanstd(ages_C, W) * fraction   # External uncertainty (est)
        σᵢ = x[i].age_sigma                 # Internal uncertainty
        σₓ = sqrt(σₑ^2 + σᵢ^2)
        x[i] = setproperty!!(x[i], :age_sigma, σₓ)
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
chronometers(data, model) = chronometers(Float64, data, model)
function chronometers(T::Type{<:AbstractFloat}, data, model)
    agesteps = floatrange(model.agesteps)
    tsteps = floatrange(model.tsteps)
    dr = haskey(model, :dr) ? model.dr : one(T)
    @assert issorted(tsteps)
    @assert tsteps == reverse(agesteps)

    haskey(data, :mineral) || @error "data must contain a column labeled `mineral`"
    mineral = data.mineral
    crystage = if haskey(data, :crystallization_age_Ma)        
        data.crystallization_age_Ma
    elseif haskey(data, :crystAge) # Legacy option
        data.crystAge
    else
        @error "data must contain a column labeled `crystallization age [Ma]`"
    end
    @assert eachindex(mineral) == eachindex(crystage)

    result = Chronometer[]
    for i in eachindex(mineral)
        first_index = 1 + round(Int,(maximum(agesteps) - crystage[i])/step(tsteps))

        if data.mineral[i] == "zircon"
            if haskey(data, :raw_He_age_Ma) && haskey(data, :raw_He_age_sigma_Ma) && !isnan(data.raw_He_age_Ma[i]/data.raw_He_age_sigma_Ma[i])
                # Modern format
                c = ZirconHe(T;
                    age = data.raw_He_age_Ma[i], 
                    age_sigma = data.raw_He_age_sigma_Ma[i], 
                    r = data.halfwidth_um[i], 
                    dr = dr, 
                    U238 = data.U238_ppm[i], 
                    Th232 = data.Th232_ppm[i], 
                    Sm147 = (haskey(data, :Sm147_ppm) && !isnan(data.Sm147_ppm[i])) ? data.Sm147_ppm[i] : 0,
                    agesteps = agesteps[first_index:end],
                )
                push!(result, c)
            elseif haskey(data, :HeAge) && haskey(data, :HeAge_sigma) && !isnan(data.HeAge[i]/data.HeAge_sigma[i])
                # Legacy format
                c = ZirconHe(T;
                    age = data.HeAge[i], 
                    age_sigma = data.HeAge_sigma[i], 
                    r = data.halfwidth[i], 
                    dr = dr, 
                    U238 = data.U[i], 
                    Th232 = data.Th[i], 
                    Sm147 = (haskey(data, :Sm) && !isnan(data.Sm[i])) ? data.Sm[i] : 0,
                    agesteps = agesteps[first_index:end],
                )
                push!(result, c)
            end

        elseif data.mineral[i] == "apatite"
            # Apatite helium
            if haskey(data, :raw_He_age_Ma) && haskey(data, :raw_He_age_sigma_Ma) && !isnan(data.raw_He_age_Ma[i]/data.raw_He_age_sigma_Ma[i])
                # Modern format
                c = ApatiteHe(T;
                    age = data.raw_He_age_Ma[i], 
                    age_sigma = data.raw_He_age_sigma_Ma[i], 
                    r = data.halfwidth_um[i], 
                    dr = dr, 
                    U238 = data.U238_ppm[i], 
                    Th232 = data.Th232_ppm[i], 
                    Sm147 = (haskey(data, :Sm147_ppm) && !isnan(data.Sm147_ppm[i])) ? data.Sm147_ppm[i] : 0,
                    agesteps = agesteps[first_index:end],
                )
                push!(result, c)
            elseif haskey(data, :HeAge) && haskey(data, :HeAge_sigma) && !isnan(data.HeAge[i]/data.HeAge_sigma[i])
                # Legacy format
                c = ApatiteHe(T;
                    age = data.HeAge[i], 
                    age_sigma = data.HeAge_sigma[i], 
                    r = data.halfwidth[i], 
                    dr = dr, 
                    U238 = data.U[i], 
                    Th232 = data.Th[i], 
                    Sm147 = (haskey(data, :Sm) && !isnan(data.Sm[i])) ? data.Sm[i] : 0,
                    agesteps = agesteps[first_index:end],
                )
                push!(result, c)
            end
            # Apatite fission track
            if haskey(data, :FT_age_Ma) && haskey(data, :FT_age_sigma_Ma) && !isnan(data.FT_age_Ma[i]/data.FT_age_sigma_Ma[i])
                c = ApatiteFT(T;
                    age = data.FT_age_Ma[i], 
                    age_sigma = data.FT_age_sigma_Ma[i], 
                    dpar = haskey(data, :dpar_um) ? data.dpar_um[i] : NaN,
                    F = haskey(data, :F_apfu) ? data.F_apfu[i] : NaN,
                    Cl = haskey(data, :Cl_apfu) ? data.Cl_apfu[i] : NaN,
                    OH = haskey(data, :OH_apfu) ? data.OH_apfu[i] : NaN,
                    rmr0 =haskey(data, :rmr0) ? data.rmr0[i] : NaN,
                    agesteps = agesteps[first_index:end],
                )
                push!(result, c)
            end
            # Apatite fission track length
            if haskey(data, :track_length_um) && haskey(data, :track_angle_degrees) && (0 < data.track_length_um[i])
                c = ApatiteTrackLength(T;
                    length = data.track_length_um[i], 
                    angle = data.track_angle_degrees[i], 
                    dpar = haskey(data, :dpar_um) ? data.dpar_um[i] : NaN,
                    F = haskey(data, :F_apfu) ? data.F_apfu[i] : NaN,
                    Cl = haskey(data, :Cl_apfu) ? data.Cl_apfu[i] : NaN,
                    OH = haskey(data, :OH_apfu) ? data.OH_apfu[i] : NaN,
                    rmr0 =haskey(data, :rmr0) ? data.rmr0[i] : NaN,
                    agesteps = agesteps[first_index:end],
                )
                push!(result, c)
            end
        end
    end

    # # Print info about samples found
    # if any(x->isa(x, ZirconHe), result)
    #     n = count(x->isa(x, ZirconHe), result)
    #     @info "found $n zircon helium samples"
    # end
    # if any(x->isa(x, ApatiteHe), result)
    #     n = count(x->isa(x, ApatiteHe), result)
    #     @info "found $n apatite helium samples"
    # end
    # if any(x->isa(x, ApatiteFT), result)
    #     n = count(x->isa(x, ApatiteFT), result)
    #     @info "found $n apatite fission track samples"
    # end
    # if any(x->isa(x, ApatiteTrackLength), result)
    #     n = count(x->isa(x, ApatiteTrackLength), result)
    #     @info "found $n apatite fission track lengths"
    # end

    isempty(result) && @error "No chronometers found"
    return unionize(result)
end

## --- model uncertainties of chronometers

function movesigma!(σcalc::AbstractVector{T}, chrons::AbstractVector{<:Chronometer}) where {T<:AbstractFloat}
    for C in (ZirconHe, ApatiteHe, ApatiteFT)
        r = abs(randn(T))
        for i in eachindex(σcalc, chrons)
            if chrons[i] isa C
                σcalc[i] *= r
            end
        end
    end
    return σcalc
end
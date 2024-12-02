# Abstract type to include any number of mineral chronometers (zircon, apatite, etc.)
abstract type Chronometer{T} end
export Chronometer

# Abstract subtypes for different categories of chronometers
abstract type FissionTrackLength{T} <: Chronometer{T} end
abstract type FissionTrackSample{T} <: Chronometer{T} end
abstract type HeliumSample{T} <: Chronometer{T} end

## --- Fission track sample types

struct ApatiteTrackLength{T<:AbstractFloat} <: FissionTrackLength{T}
    l::T                    # [um]
    theta::T                # [degrees]
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
export ApatiteTrackLength

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
export ApatiteFT

## --- Helium sample types
"""
```julia
ZirconHe(r, dr, U238, Th232, [Sm147], agesteps::AbstractVector)
```
Construct a `ZirconHe` object
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
    r238U = U238.*ones(T, size(rsteps)) # [PPMw]
    r235U = T.(r238U/137.818) # [PPMw]
    r232Th = Th232 .* ones(T, size(rsteps)) # [PPMw]
    r147Sm = Sm147 .* ones(T, size(rsteps)) # [PPMw]

    # Convert to atoms per gram
    r238U *= 6.022E23 / 1E6 / 238
    r235U *= 6.022E23 / 1E6 / 235
    r232Th *= 6.022E23 / 1E6 / 232
    r147Sm *= 6.022E23 / 1E6 / 147

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance
    dint = zeros(T, length(redges) - 1)
    
    r238UHe = zeros(T, size(r238U))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii238U)
            # Effective radial alpha deposition from 238U
            intersectiondensity!(dint,redges,relvolumes,alpharadii238U[i],rsteps[ri])
            @. r238UHe += relvolumes[ri] * dint * r238U[ri]
        end
    end

    # Effective radial alpha deposition from U-235
    r235UHe = zeros(T, size(r235U))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii235U)
            # Effective radial alpha deposition from 235U
            intersectiondensity!(dint, redges,relvolumes,alpharadii235U[i],rsteps[ri])
            @. r235UHe += relvolumes[ri] * dint * r235U[ri]
        end
    end

    # Effective radial alpha deposition from Th-232
    r232ThHe = zeros(T, size(r232Th))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii232Th)
            # Effective radial alpha deposition from 232Th
            intersectiondensity!(dint, redges,relvolumes,alpharadii232Th[i],rsteps[ri])
            @. r232ThHe += relvolumes[ri] * dint * r232Th[ri]
        end
    end

    # Effective radial alpha deposition from Sm-147
    r147SmHe = zeros(T, size(r147Sm))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii147Sm)
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
export ZirconHe


"""
```julia
ApatiteHe(r, dr, U238, Th232, [Sm147], agesteps::AbstractVector)
```
Construct an `ApatiteHe` object
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
    r232Th = fill(T(Th232), size(rsteps))    # [PPMw]
    r147Sm = fill(T(Sm147), size(rsteps))    # [PPMw]

    # Convert to atoms per gram
    r238U *= 6.022E23 / 1E6 / 238
    r235U *= 6.022E23 / 1E6 / 235
    r232Th *= 6.022E23 / 1E6 / 232
    r147Sm *= 6.022E23 / 1E6 / 147

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance
    dint = zeros(T, length(redges) - 1)

    # Effective radial alpha deposition from U-238
    r238UHe = zeros(T, size(r238U))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii238U)
            # Effective radial alpha deposition from 238U
            intersectiondensity!(dint,redges,relvolumes,alpharadii238U[i],rsteps[ri])
            @. r238UHe += relvolumes[ri] * dint * r238U[ri]
        end
    end

    # Effective radial alpha deposition from U-235
    r235UHe = zeros(T, size(r235U))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii235U)
            # Effective radial alpha deposition from 235U
            intersectiondensity!(dint, redges,relvolumes,alpharadii235U[i],rsteps[ri])
            @. r235UHe += relvolumes[ri] * dint * r235U[ri]
        end
    end

    # Effective radial alpha deposition from Th-232
    r232ThHe = zeros(T, size(r232Th))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii232Th)
            # Effective radial alpha deposition from 232Th
            intersectiondensity!(dint, redges,relvolumes,alpharadii232Th[i],rsteps[ri])
            @. r232ThHe += relvolumes[ri] * dint * r232Th[ri]
        end
    end

    # Effective radial alpha deposition from Sm-147
    r147SmHe = zeros(T, size(r147Sm))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii147Sm)
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
export ApatiteHe


"""
```julia
chronometer([T=Float64], data, model)
```
Construct a vector of `Chronometer` objects given a dataset `data`
and model parameters `model`
"""
chronometer(data, model) = chronometer(Float64, data, model)
function chronometer(T::Type{<:AbstractFloat}, data, model)
    agesteps = floatrange(model.agesteps)
    tsteps = floatrange(model.tsteps)
    dr = haskey(model, :dr) ? model.dr : one(T)
    @assert issorted(tsteps)
    @assert tsteps == reverse(agesteps)

    haskey(data, :mineral) || @error "data must contain a column labeled `mineral`"
    mineral = data.mineral
    haskey(data, :crystallization_age_Ma) || @error "data must contain a column labeled `crystallization age [Ma]`"
    crystage = data.crystallization_age_Ma
    @assert eachindex(data.mineral) == eachindex(data.crystallization_age_Ma)

    result = Chronometer[]
    for i in eachindex(mineral)
        first_index = 1 + round(Int,(maximum(agesteps) - crystage[i])/step(tsteps))

        if data.mineral[i] == "zircon"
            if haskey(data, :HeAge) && haskey(data, :HeAge_sigma) && !isnan(data.HeAge[i]/data.HeAge_sigma[i])
                c = ZirconHe(T;
                    age = data.raw_He_age_Ma[i], 
                    age_sigma = data.raw_He_age_sigma_Ma[i], 
                    r = data.halfwidth_um[i], 
                    dr = dr, 
                    U238 = data.U238_ppm[i], 
                    Th232 = data.Th232_ppm[i], 
                    Sm147 = (haskey(data, :Sm147_ppm) && isnan(data.Sm147_ppm[i])) ? data.Sm147_ppm[i] : 0,
                    agesteps = agesteps[first_index:end],
                )
                push!(result, c)
            end

        elseif data.mineral[i] == "apatite"
            if haskey(data, :raw_He_age_Ma) && haskey(data, :raw_He_age_sigma_Ma) && !isnan(data.raw_He_age_Ma[i]/data.raw_He_age_sigma_Ma[i])
                c = ApatiteHe(T;
                    age = data.raw_He_age_Ma[i], 
                    age_sigma = data.raw_He_age_sigma_Ma[i], 
                    r = data.halfwidth_um[i], 
                    dr = dr, 
                    U238 = data.U238_ppm[i], 
                    Th232 = data.Th232_ppm[i], 
                    Sm147 = (haskey(data, :Sm147_ppm) && isnan(data.Sm147_ppm[i])) ? data.Sm147_ppm[i] : 0,
                    agesteps = agesteps[first_index:end],
                )
                push!(result, c)
            end
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
    return unionize(result)
end
export chronometer
# Abstract type to include any number of mineral chronometers (zircon, apatite, etc.)
abstract type Mineral{T} end


"""
```julia
ZirconHe(r, dr, Uppm, Th232ppm, [Sm147ppm], dt, agesteps::AbstractVector)
```
Construct a `ZirconHe` object
"""
# Concretely-typed immutable struct to hold information about a single zircon (Helium) crystal
struct ZirconHe{T<:Number} <: Mineral{T}
    dt::T
    agesteps::Vector{T}
    tsteps::Vector{T}
    ntsteps::Int
    dr::T
    rsteps::Vector{T}
    redges::Vector{T}
    relvolumes::Vector{T}
    nrsteps::Int
    r238U::Vector{T}
    r235U::Vector{T}
    r232Th::Vector{T}
    r147Sm::Vector{T}
    r238UHe::Vector{T}
    r235UHe::Vector{T}
    r232ThHe::Vector{T}
    r147SmHe::Vector{T}
    alphadeposition::Matrix{T}
    alphadamage::Matrix{T}
    annealeddamage::Matrix{T}
    u::Matrix{T}
    β::Vector{T}
    Dz::Vector{T}
    DN17::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{Float64, Tridiagonal{Float64, Vector{Float64}}, Vector{Int64}}
    y::Vector{T}
end
# Constructor for the ZirconHe type, given grain radius, U and Th concentrations and t-T discretization information
ZirconHe(r::T, dr::Number, Uppm::T, Th232ppm::T, dt::Number, agesteps::AbstractVector{T}) where T<:Number = ZirconHe(r, dr, Uppm, Th232ppm, zero(T), dt, agesteps)
function ZirconHe(r::T, dr::Number, Uppm::T, Th232ppm::T, Sm147ppm::T, dt::Number, agesteps::AbstractVector{T}) where T<:Number
    
    # Temporal discretization
    tsteps = reverse(agesteps)
    ntsteps = length(tsteps) # Number of time steps

    # crystal size and spatial discretization
    rsteps = Array{T}(0+dr/2 : dr: r-dr/2)
    redges = Array{T}(0 : dr : r) # Edges of each radius element
    nrsteps = length(rsteps)+2 # number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    relvolumes = (redges[2:end].^3 - redges[1:end-1].^3)/redges[end]^3 # Relative volume fraction of spherical shell corresponding to each radius element

    # Alpha stopping distances for each isotope in each decay chain, from
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alpharadii238U = (11.78, 14.09, 13.73, 14.13, 17.32, 16.69, 28.56, 16.48,)
    alpharadii235U = (12.58, 15.04, 19.36, 18.06, 23.07, 26.87, 22.47,)
    alpharadii232Th = (10.99, 16.67, 18.16, 17.32, 23.61, 29.19,)
    # Ketchem et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alpharadii147Sm = (4.76,)

    # Decay constants
    λ235U = log(2)/(7.0381*10^8)*10^6 # [1/Myr] Jaffey et al. (1971)
    λ238U = log(2)/(4.4683*10^9)*10^6 # [1/Myr] Jaffey et al. (1971)
    λ232Th = log(2)/(1.405*10^10)*10^6 # [1/Myr]
    λ147Sm = log(2)/(1.07*10^11)*10^6 # [1/Myr]

    # Observed radial HPE profiles at present day
    r238U = Uppm.*ones(T, size(rsteps)) # [PPMw]
    r235U = T.(r238U/137.818) # [PPMw]
    r232Th = Th232ppm .* ones(T, size(rsteps)) # [PPMw]
    r147Sm = Sm147ppm .* ones(T, size(rsteps)) # [PPMw]

    # Convert to atoms per gram
    r238U *= 6.022E23 / 1E6 / 238
    r235U *= 6.022E23 / 1E6 / 235
    r232Th *= 6.022E23 / 1E6 / 232
    r147Sm *= 6.022E23 / 1E6 / 147

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance
    dint = zeros(T, length(redges) - 1)

    
    r238UHe = zeros(size(r238U))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii238U)
            # Effective radial alpha deposition from 238U
            intersectiondensity!(dint,redges,relvolumes,alpharadii238U[i],rsteps[ri])
            @turbo @. r238UHe += relvolumes[ri] * dint * r238U[ri]
        end
    end

    # Effective radial alpha deposition from U-235
    r235UHe = zeros(size(r235U))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii235U)
            # Effective radial alpha deposition from 235U
            intersectiondensity!(dint, redges,relvolumes,alpharadii235U[i],rsteps[ri])
            @turbo @. r235UHe += relvolumes[ri] * dint * r235U[ri]
        end
    end

    # Effective radial alpha deposition from Th-232
    r232ThHe = zeros(size(r232Th))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii232Th)
            # Effective radial alpha deposition from 232Th
            intersectiondensity!(dint, redges,relvolumes,alpharadii232Th[i],rsteps[ri])
            @turbo @. r232ThHe += relvolumes[ri] * dint * r232Th[ri]
        end
    end

    # Effective radial alpha deposition from Sm-147
    r147SmHe = zeros(size(r147Sm))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii147Sm)
            intersectiondensity!(dint, redges,relvolumes,alpharadii147Sm[i],rsteps[ri])
            @turbo @. r147SmHe += relvolumes[ri] * dint * r147Sm[ri]
        end
    end    

    # Alpha decay recoil damage
    r238Udam = 8*r238U # No smoothing of alpha damage, 8 alphas per 238U
    r235Udam = 7*r235U # No smoothing of alpha damage, 7 alphas per 235U
    r232Thdam = 6*r232Th # No smoothing of alpha damage, 6 alphas per 232 Th
    r147Smdam = 1*r147Sm # No smoothing of alpha damage, 1 alpha per 147 Sm

    # Calculate corrected alpha deposition and recoil damage each time step for each radius
    dt_2 = dt/2
    decay = Array{T}(undef, ntsteps)
    buffer = zeros(T, ntsteps, nrsteps-2)
    # Allocate deposition and damage arrays
    alphadeposition = zeros(T, ntsteps, nrsteps-2)
    alphadamage = zeros(T, ntsteps, nrsteps-2)
    # U-238
    @turbo @. decay = exp(λ238U*(agesteps + dt_2)) - exp(λ238U*(agesteps - dt_2))
    mul!(buffer, decay, r238UHe')
    @turbo @. alphadeposition += buffer
    mul!(buffer, decay, r238Udam')
    @turbo @. alphadamage += buffer
    # U-232
    @turbo @. decay = exp(λ235U*(agesteps + dt_2)) - exp(λ235U*(agesteps - dt_2))
    mul!(buffer, decay, r235UHe')
    @turbo @. alphadeposition += buffer
    mul!(buffer, decay, r235Udam')
    @turbo @. alphadamage += buffer
    # Th-232
    @turbo @. decay = exp(λ232Th*(agesteps + dt_2)) - exp(λ232Th*(agesteps - dt_2))
    mul!(buffer, decay, r232ThHe')
    @turbo @. alphadeposition += buffer
    mul!(buffer, decay, r232Thdam')
    @turbo @. alphadamage += buffer
    # Sm-147
    @turbo @. decay = exp(λ147Sm*(agesteps + dt_2)) - exp(λ147Sm*(agesteps - dt_2))
    mul!(buffer, decay, r147SmHe')
    @turbo @. alphadeposition += buffer
    mul!(buffer, decay, r147Smdam')
    @turbo @. alphadamage += buffer

    # Allocate additional variables that will be needed for Crank-Nicholson
    annealeddamage = similar(alphadamage)
    β = Array{T}(undef, nrsteps) # First row of annealeddamage

    # Allocate arrays for diffusivities
    Dz = Array{T}(undef, ntsteps)
    DN17 = Array{T}(undef, ntsteps)

    # Allocate output matrix for all timesteps
    u = Array{T}(undef, nrsteps, ntsteps)

    # Allocate variables for tridiagonal matrix and RHS
    dl = ones(T, nrsteps-1)    # Sub-diagonal row
    d = ones(T, nrsteps)       # Diagonal
    du = ones(T, nrsteps-1)    # Supra-diagonal row
    du2 = ones(T, nrsteps-2)   # sup-sup-diagonal row for pivoting

    # Tridiagonal matrix for LHS of Crank-Nicholson equation with regular grid cells
    A = Tridiagonal(dl, d, du, du2)
    F = lu(A, allowsingular=true)

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = Array{T}(undef, nrsteps)

    return ZirconHe{T}(
        dt,
        agesteps,
        tsteps,
        ntsteps,
        dr,
        rsteps,
        redges,
        relvolumes,
        nrsteps,
        r238U,
        r235U,
        r232Th,
        r147Sm,
        r238UHe,
        r235UHe,
        r232ThHe,
        r147SmHe,
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
ApatiteHe(r, dr, Uppm, Th232ppm, [Sm147ppm], dt, agesteps::AbstractVector)
```
Construct an `ApatiteHe` object
"""
# Concretely-typed immutable struct to hold information about a single apatite (Helium) crystal
struct ApatiteHe{T<:Number} <: Mineral{T}
    dt::T
    agesteps::Vector{T}
    tsteps::Vector{T}
    ntsteps::Int
    dr::T
    rsteps::Vector{T}
    redges::Vector{T}
    relvolumes::Vector{T}
    nrsteps::Int
    r238U::Vector{T}
    r235U::Vector{T}
    r232Th::Vector{T}
    r147Sm::Vector{T}
    r238UHe::Vector{T}
    r235UHe::Vector{T}
    r232ThHe::Vector{T}
    r147SmHe::Vector{T}
    alphadeposition::Matrix{T}
    alphadamage::Matrix{T}
    annealeddamage::Matrix{T}
    u::Matrix{T}
    β::Vector{T}
    DL::Vector{T}
    Dtrap::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    F::LU{Float64, Tridiagonal{Float64, Vector{Float64}}, Vector{Int64}}
    y::Vector{T}
end
# Constructor for the ApatiteHe type, given grain radius, U and Th concentrations and t-T discretization information
ApatiteHe(r::T, dr::Number, Uppm::T, Th232ppm::T, dt::Number, agesteps::AbstractVector{T}) where T<:Number = ApatiteHe(r, dr, Uppm, Th232ppm, zero(T), dt, agesteps)
function ApatiteHe(r::T, dr::Number, Uppm::T, Th232ppm::T, Sm147ppm::T, dt::Number, agesteps::AbstractVector{T}) where T<:Number

    # Temporal discretization
    tsteps = reverse(agesteps)
    ntsteps = length(tsteps) # Number of time steps

    # crystal size and spatial discretization
    rsteps = Array{T}(0+dr/2 : dr: r-dr/2)
    redges = Array{T}(0 : dr : r) # Edges of each radius element
    nrsteps = length(rsteps)+2 # number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    relvolumes = (redges[2:end].^3 - redges[1:end-1].^3)/redges[end]^3 # Relative volume fraction of spherical shell corresponding to each radius element

    # Alpha stopping distances for each isotope in each decay chain, from
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alpharadii238U = (13.54, 16.26, 15.84, 16.31, 20.09, 22.89, 33.39, 19.10,)
    alpharadii235U = (14.48, 17.39, 22.5, 20.97, 26.89, 31.40, 26.18,)
    alpharadii232Th = (12.60, 19.32, 21.08, 20.09, 27.53, 34.14,)
    # Ketchem et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alpharadii147Sm = (5.93,)

    # Jaffey decay constants
    λ235U = log(2)/(7.0381*10^8)*10^6 # [1/Myr]
    λ238U = log(2)/(4.4683*10^9)*10^6 # [1/Myr]
    λ232Th = log(2)/(1.405*10^10)*10^6 # [1/Myr]
    λ147Sm = log(2)/(1.07*10^11)*10^6 # [1/Myr]

    # Observed radial HPE profiles at present day
    r238U = Uppm.*ones(T, size(rsteps)) # PPM
    r235U = T.(r238U/137.818) #PPM
    r232Th = Th232ppm .* ones(T, size(rsteps)) #PPM
    r147Sm = Sm147ppm .* ones(T, size(rsteps)) # [PPMw]

    # Convert to atoms per gram
    r238U *= 6.022E23 / 1E6 / 238
    r235U *= 6.022E23 / 1E6 / 235
    r232Th *= 6.022E23 / 1E6 / 232
    r147Sm *= 6.022E23 / 1E6 / 147

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance
    dint = zeros(T, length(redges) - 1)

    # Effective radial alpha deposition from U-238
    r238UHe = zeros(size(r238U))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii238U)
            # Effective radial alpha deposition from 238U
            intersectiondensity!(dint,redges,relvolumes,alpharadii238U[i],rsteps[ri])
            @turbo @. r238UHe += relvolumes[ri] * dint * r238U[ri]
        end
    end

    # Effective radial alpha deposition from U-235
    r235UHe = zeros(size(r235U))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii235U)
            # Effective radial alpha deposition from 235U
            intersectiondensity!(dint, redges,relvolumes,alpharadii235U[i],rsteps[ri])
            @turbo @. r235UHe += relvolumes[ri] * dint * r235U[ri]
        end
    end

    # Effective radial alpha deposition from Th-232
    r232ThHe = zeros(size(r232Th))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii232Th)
            # Effective radial alpha deposition from 232Th
            intersectiondensity!(dint, redges,relvolumes,alpharadii232Th[i],rsteps[ri])
            @turbo @. r232ThHe += relvolumes[ri] * dint * r232Th[ri]
        end
    end

    # Effective radial alpha deposition from Sm-147
    r147SmHe = zeros(size(r147Sm))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alpharadii147Sm)
            intersectiondensity!(dint, redges,relvolumes,alpharadii147Sm[i],rsteps[ri])
            @turbo @. r147SmHe += relvolumes[ri] * dint * r147Sm[ri]
        end
    end

    # Alpha decay recoil damage
    r238Udam = 8*r238U # No smoothing of alpha damage, 8 alphas per 238U
    r235Udam = 7*r235U # No smoothing of alpha damage, 7 alphas per 235U
    r232Thdam = 6*r232Th # No smoothing of alpha damage, 6 alphas per 232 Th
    r147Smdam = 1*r147Sm # No smoothing of alpha damage, 1 alpha per 147 Sm

    # Calculate corrected alpha deposition and recoil damage each time step for each radius
    dt_2 = dt/2
    decay = Array{T}(undef, ntsteps)
    buffer = zeros(T, ntsteps, nrsteps-2)
    # Allocate deposition and damage arrays
    alphadeposition = zeros(T, ntsteps, nrsteps-2)
    alphadamage = zeros(T, ntsteps, nrsteps-2)
    # U-238
    @turbo @. decay = exp(λ238U*(agesteps + dt_2)) - exp(λ238U*(agesteps - dt_2))
    mul!(buffer, decay, r238UHe')
    @turbo @. alphadeposition += buffer
    mul!(buffer, decay, r238Udam')
    @turbo @. alphadamage += buffer
    # U-232
    @turbo @. decay = exp(λ235U*(agesteps + dt_2)) - exp(λ235U*(agesteps - dt_2))
    mul!(buffer, decay, r235UHe')
    @turbo @. alphadeposition += buffer
    mul!(buffer, decay, r235Udam')
    @turbo @. alphadamage += buffer
    # Th-232
    @turbo @. decay = exp(λ232Th*(agesteps + dt_2)) - exp(λ232Th*(agesteps - dt_2))
    mul!(buffer, decay, r232ThHe')
    @turbo @. alphadeposition += buffer
    mul!(buffer, decay, r232Thdam')
    @turbo @. alphadamage += buffer
    # Sm-147
    @turbo @. decay = exp(λ147Sm*(agesteps + dt_2)) - exp(λ147Sm*(agesteps - dt_2))
    mul!(buffer, decay, r147SmHe')
    @turbo @. alphadeposition += buffer
    mul!(buffer, decay, r147Smdam')
    @turbo @. alphadamage += buffer

    # Allocate additional variables that will be needed for Crank-Nicholson
    annealeddamage = similar(alphadamage)
    β = Array{T}(undef, nrsteps) # First row of annealeddamage

    # Allocate arrays for diffusivities
    DL = Array{T}(undef, ntsteps)
    Dtrap = Array{T}(undef, ntsteps)

    # Allocate output matrix for all timesteps
    u = Array{T}(undef, nrsteps, ntsteps)

    # Allocate variables for tridiagonal matrix and RHS
    dl = ones(T, nrsteps-1)    # Sub-diagonal row
    d = ones(T, nrsteps)       # Diagonal
    du = ones(T, nrsteps-1)    # Supra-diagonal row
    du2 = ones(T, nrsteps-2)   # sup-sup-diagonal row for pivoting

    # Tridiagonal matrix for LHS of Crank-Nicholson equation with regular grid cells
    A = Tridiagonal(dl, d, du, du2)
    F = lu(A, allowsingular=true)

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = Array{T}(undef, nrsteps)

    return ApatiteHe{T}(
        dt,
        agesteps,
        tsteps,
        ntsteps,
        dr,
        rsteps,
        redges,
        relvolumes,
        nrsteps,
        r238U,
        r235U,
        r232Th,
        r147Sm,
        r238UHe,
        r235UHe,
        r232ThHe,
        r147SmHe,
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

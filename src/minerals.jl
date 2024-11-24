"""
```julia
intersectionfraction(r₁, r₂, d)
```
Calculate the fraction of the surface area of a sphere s2 with radius `r₂`
that intersects the interior of a sphere s1 of radius `r₁` if the two are
separated by distance `d`. Assumptions: `r₁`>0, `r₂`>0, `d`>=0
"""
function intersectionfraction(r₁::T, r₂::T, d::T) where T <: AbstractFloat
    dmax = r₁+r₂
    dmin = abs(r₁-r₂)
    if d > dmax ## If separated by more than dmax, there is no intersection
        omega = zero(T)
    elseif d > dmin
        # X is the radial distance between the center of s2 and the interseciton plane
        # See e.g., http://mathworld.wolfram.com/Sphere-SphereIntersection.html
        # x = (d^2 - r₁^2 + r₂^2) / (2 * d)

        # Let omega be is the solid angle of intersection normalized by 4pi,
        # where the solid angle of a cone is 2pi*(1-cos(theta)) and cos(theta)
        # is adjacent/hypotenuse = x/r₂.
        # omega = (1 - x/r₂)/2

        # Rearranged for optimization:
        omega = one(T)/2 - (d^2 - r₁^2 + r₂^2) / (4 * d * r₂)

    elseif r₁<r₂ # If r₁ is entirely within r₂
        omega = zero(T)
    else
        omega = one(T) # If r₂ is entirely within r₁
    end
    return omega
end
export intersectionfraction


"""
```julia
intersectiondensity(rEdges::Vector, rVolumes::Vector, ralpha, d)
```
Calculate the volume-nomalized fractional intersection density of an alpha
stopping sphere of radius `ralpha` with each concentric shell (with shell edges
`rEdges` and relative volumes `rVolumes`) of a spherical crystal where the
two are separated by distance `d`
"""
function intersectiondensity(rEdges::Vector{T}, rVolumes::Vector{T}, ralpha::T, d::T) where T <: AbstractFloat
    dInt = Array{T}(undef, length(rEdges) - 1)
    intersectiondensity!(dInt, rEdges, rVolumes, ralpha, d)
end
function intersectiondensity!(dInt::Vector{T}, rEdges::Vector{T}, rVolumes::Vector{T}, ralpha::T, d::T) where T <: AbstractFloat
    n = length(rEdges) - 1
    fIntLast = intersectionfraction(first(rEdges),ralpha,d)
    @inbounds for i ∈ 1:n
        # Integrated intersection fraction for each concentric sphere (rEdges) of crystal
        fInt = intersectionfraction(rEdges[i+1],ralpha,d)

        # Intersection fraction for each spherical shell of the crystal (subtracting
        # one concentric sphere from the next) normalized by shell volume
        dInt[i] = (fInt-fIntLast) / rVolumes[i]
        fIntLast = fInt
    end
    return dInt
end
export intersectiondensity

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
    rEdges::Vector{T}
    relVolumes::Vector{T}
    nrsteps::Int
    r238U::Vector{T}
    r235U::Vector{T}
    r232Th::Vector{T}
    r147Sm::Vector{T}
    r238UHe::Vector{T}
    r235UHe::Vector{T}
    r232ThHe::Vector{T}
    r147SmHe::Vector{T}
    alphaDeposition::Matrix{T}
    alphaDamage::Matrix{T}
    annealedDamage::Matrix{T}
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
    rEdges = Array{T}(0 : dr : r) # Edges of each radius element
    nrsteps = length(rsteps)+2 # number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    relVolumes = (rEdges[2:end].^3 - rEdges[1:end-1].^3)/rEdges[end]^3 # Relative volume fraction of spherical shell corresponding to each radius element

    # Alpha stopping distances for each isotope in each decay chain, from
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alphaRadii238U = (11.78, 14.09, 13.73, 14.13, 17.32, 16.69, 28.56, 16.48,)
    alphaRadii235U = (12.58, 15.04, 19.36, 18.06, 23.07, 26.87, 22.47,)
    alphaRadii232Th = (10.99, 16.67, 18.16, 17.32, 23.61, 29.19,)
    # Ketchem et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alphaRadii147Sm = (4.76,)

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
    dInt = zeros(T, length(rEdges) - 1)

    
    r238UHe = zeros(size(r238U))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alphaRadii238U)
            # Effective radial alpha deposition from 238U
            intersectiondensity!(dInt,rEdges,relVolumes,alphaRadii238U[i],rsteps[ri])
            @turbo @. r238UHe += relVolumes[ri] * dInt * r238U[ri]
        end
    end

    # Effective radial alpha deposition from U-235
    r235UHe = zeros(size(r235U))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alphaRadii235U)
            # Effective radial alpha deposition from 235U
            intersectiondensity!(dInt, rEdges,relVolumes,alphaRadii235U[i],rsteps[ri])
            @turbo @. r235UHe += relVolumes[ri] * dInt * r235U[ri]
        end
    end

    # Effective radial alpha deposition from Th-232
    r232ThHe = zeros(size(r232Th))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alphaRadii232Th)
            # Effective radial alpha deposition from 232Th
            intersectiondensity!(dInt, rEdges,relVolumes,alphaRadii232Th[i],rsteps[ri])
            @turbo @. r232ThHe += relVolumes[ri] * dInt * r232Th[ri]
        end
    end

    # Effective radial alpha deposition from Sm-147
    r147SmHe = zeros(size(r147Sm))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alphaRadii147Sm)
            intersectiondensity!(dInt, rEdges,relVolumes,alphaRadii147Sm[i],rsteps[ri])
            @turbo @. r147SmHe += relVolumes[ri] * dInt * r147Sm[ri]
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
    alphaDeposition = zeros(T, ntsteps, nrsteps-2)
    alphaDamage = zeros(T, ntsteps, nrsteps-2)
    # U-238
    @turbo @. decay = exp(λ238U*(agesteps + dt_2)) - exp(λ238U*(agesteps - dt_2))
    mul!(buffer, decay, r238UHe')
    @turbo @. alphaDeposition += buffer
    mul!(buffer, decay, r238Udam')
    @turbo @. alphaDamage += buffer
    # U-232
    @turbo @. decay = exp(λ235U*(agesteps + dt_2)) - exp(λ235U*(agesteps - dt_2))
    mul!(buffer, decay, r235UHe')
    @turbo @. alphaDeposition += buffer
    mul!(buffer, decay, r235Udam')
    @turbo @. alphaDamage += buffer
    # Th-232
    @turbo @. decay = exp(λ232Th*(agesteps + dt_2)) - exp(λ232Th*(agesteps - dt_2))
    mul!(buffer, decay, r232ThHe')
    @turbo @. alphaDeposition += buffer
    mul!(buffer, decay, r232Thdam')
    @turbo @. alphaDamage += buffer
    # Sm-147
    @turbo @. decay = exp(λ147Sm*(agesteps + dt_2)) - exp(λ147Sm*(agesteps - dt_2))
    mul!(buffer, decay, r147SmHe')
    @turbo @. alphaDeposition += buffer
    mul!(buffer, decay, r147Smdam')
    @turbo @. alphaDamage += buffer

    # Allocate additional variables that will be needed for Crank-Nicholson
    annealedDamage = similar(alphaDamage)
    β = Array{T}(undef, nrsteps) # First row of annealedDamage

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
        rEdges,
        relVolumes,
        nrsteps,
        r238U,
        r235U,
        r232Th,
        r147Sm,
        r238UHe,
        r235UHe,
        r232ThHe,
        r147SmHe,
        alphaDeposition,
        alphaDamage,
        annealedDamage,
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
    rEdges::Vector{T}
    relVolumes::Vector{T}
    nrsteps::Int
    r238U::Vector{T}
    r235U::Vector{T}
    r232Th::Vector{T}
    r147Sm::Vector{T}
    r238UHe::Vector{T}
    r235UHe::Vector{T}
    r232ThHe::Vector{T}
    r147SmHe::Vector{T}
    alphaDeposition::Matrix{T}
    alphaDamage::Matrix{T}
    annealedDamage::Matrix{T}
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
    rEdges = Array{T}(0 : dr : r) # Edges of each radius element
    nrsteps = length(rsteps)+2 # number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    relVolumes = (rEdges[2:end].^3 - rEdges[1:end-1].^3)/rEdges[end]^3 # Relative volume fraction of spherical shell corresponding to each radius element

    # Alpha stopping distances for each isotope in each decay chain, from
    # Farley et al. (1996), doi: 10.1016/S0016-7037(96)00193-7
    alphaRadii238U = (13.54, 16.26, 15.84, 16.31, 20.09, 22.89, 33.39, 19.10,)
    alphaRadii235U = (14.48, 17.39, 22.5, 20.97, 26.89, 31.40, 26.18,)
    alphaRadii232Th = (12.60, 19.32, 21.08, 20.09, 27.53, 34.14,)
    # Ketchem et al. (2011), doi: 10.1016/j.gca.2011.10.011
    alphaRadii147Sm = (5.93,)

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
    dInt = zeros(T, length(rEdges) - 1)

    # Effective radial alpha deposition from U-238
    r238UHe = zeros(size(r238U))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alphaRadii238U)
            # Effective radial alpha deposition from 238U
            intersectiondensity!(dInt,rEdges,relVolumes,alphaRadii238U[i],rsteps[ri])
            @turbo @. r238UHe += relVolumes[ri] * dInt * r238U[ri]
        end
    end

    # Effective radial alpha deposition from U-235
    r235UHe = zeros(size(r235U))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alphaRadii235U)
            # Effective radial alpha deposition from 235U
            intersectiondensity!(dInt, rEdges,relVolumes,alphaRadii235U[i],rsteps[ri])
            @turbo @. r235UHe += relVolumes[ri] * dInt * r235U[ri]
        end
    end

    # Effective radial alpha deposition from Th-232
    r232ThHe = zeros(size(r232Th))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alphaRadii232Th)
            # Effective radial alpha deposition from 232Th
            intersectiondensity!(dInt, rEdges,relVolumes,alphaRadii232Th[i],rsteps[ri])
            @turbo @. r232ThHe += relVolumes[ri] * dInt * r232Th[ri]
        end
    end

    # Effective radial alpha deposition from Sm-147
    r147SmHe = zeros(size(r147Sm))
    @inbounds for ri = 1:length(rsteps)
        for i=1:length(alphaRadii147Sm)
            intersectiondensity!(dInt, rEdges,relVolumes,alphaRadii147Sm[i],rsteps[ri])
            @turbo @. r147SmHe += relVolumes[ri] * dInt * r147Sm[ri]
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
    alphaDeposition = zeros(T, ntsteps, nrsteps-2)
    alphaDamage = zeros(T, ntsteps, nrsteps-2)
    # U-238
    @turbo @. decay = exp(λ238U*(agesteps + dt_2)) - exp(λ238U*(agesteps - dt_2))
    mul!(buffer, decay, r238UHe')
    @turbo @. alphaDeposition += buffer
    mul!(buffer, decay, r238Udam')
    @turbo @. alphaDamage += buffer
    # U-232
    @turbo @. decay = exp(λ235U*(agesteps + dt_2)) - exp(λ235U*(agesteps - dt_2))
    mul!(buffer, decay, r235UHe')
    @turbo @. alphaDeposition += buffer
    mul!(buffer, decay, r235Udam')
    @turbo @. alphaDamage += buffer
    # Th-232
    @turbo @. decay = exp(λ232Th*(agesteps + dt_2)) - exp(λ232Th*(agesteps - dt_2))
    mul!(buffer, decay, r232ThHe')
    @turbo @. alphaDeposition += buffer
    mul!(buffer, decay, r232Thdam')
    @turbo @. alphaDamage += buffer
    # Sm-147
    @turbo @. decay = exp(λ147Sm*(agesteps + dt_2)) - exp(λ147Sm*(agesteps - dt_2))
    mul!(buffer, decay, r147SmHe')
    @turbo @. alphaDeposition += buffer
    mul!(buffer, decay, r147Smdam')
    @turbo @. alphaDamage += buffer

    # Allocate additional variables that will be needed for Crank-Nicholson
    annealedDamage = similar(alphaDamage)
    β = Array{T}(undef, nrsteps) # First row of annealedDamage

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
        rEdges,
        relVolumes,
        nrsteps,
        r238U,
        r235U,
        r232Th,
        r147Sm,
        r238UHe,
        r235UHe,
        r232ThHe,
        r147SmHe,
        alphaDeposition,
        alphaDamage,
        annealedDamage,
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

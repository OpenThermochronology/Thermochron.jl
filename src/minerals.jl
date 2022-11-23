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

abstract type Mineral end

struct Zircon{T<:Number} <: Mineral
    dt::T
    ageSteps::Vector{T}
    tSteps::Vector{T}
    ntSteps::Int
    dr::T
    rSteps::Vector{T}
    rEdges::Vector{T}
    relVolumes::Vector{T}
    nrSteps::Int
    r238U::Vector{T}
    r235U::Vector{T}
    r232Th::Vector{T}
    r238UHe::Vector{T}
    r235UHe::Vector{T}
    r232ThHe::Vector{T}
    alphaDeposition::Matrix{T}
    alphaDamage::Matrix{T}
    annealedDamage::Matrix{T}
    u::Matrix{T}
    β::Vector{T}
    Dz::Vector{T}
    DN17::Vector{T}
    A::Tridiagonal{T, Vector{T}}
    ipiv::Vector{LinearAlgebra.BlasInt}
    y::Vector{T}
end

function Zircon(r::T, dr::Number, Uppm::T, Thppm::T, dt::Number, ageSteps::AbstractVector{T}) where T<:Number

    # Temporal discretization
    tSteps = reverse(ageSteps)
    ntSteps = length(tSteps) # Number of time steps

    # Crystal size and spatial discretization
    rSteps = Array{T}(0+dr/2 : dr: r-dr/2)
    rEdges = Array{T}(0 : dr : r) # Edges of each radius element
    nrSteps = length(rSteps)+2 # number of radial grid points -- note 2 implict points: one at negative radius, one outside grain
    relVolumes = (rEdges[2:end].^3 - rEdges[1:end-1].^3)/rEdges[end]^3 # Relative volume fraction of spherical shell corresponding to each radius element

    # Alpha stopping distances for each isotope in each decay chain, from
    # Farley et al., 1996
    alphaRadii238U = (11.78, 14.09, 13.73, 14.13, 17.32, 16.69, 28.56, 16.48,)
    alphaRadii235U = (12.58, 15.04, 19.36, 18.06, 23.07, 26.87, 22.47,)
    alphaRadii232Th = (10.99, 16.67, 18.16, 17.32, 23.61, 29.19,)

    # Jaffey decay constants
    λ235U = log(2)/(7.0381*10^8)*10^6 # [1/Myr]
    λ238U = log(2)/(4.4683*10^9)*10^6 # [1/Myr]
    λ232Th = log(2)/(1.405*10^10)*10^6 # [1/Myr]

    # Observed radial HPE profiles at present day
    r238U = Uppm.*ones(T, size(rSteps)) # PPM
    r235U = T.(r238U/137.818) #PPM
    r232Th = Thppm .* ones(T, size(rSteps)) #PPM

    # Convert to atoms per gram
    r238U *= 6.022E23 / 1E6 / 238
    r235U *= 6.022E23 / 1E6 / 235
    r232Th *= 6.022E23 / 1E6 / 232

    # Calculate effective He deposition for each decay chain, corrected for alpha
    # stopping distance
    dInt = Array{Float64}(undef, length(rEdges) - 1)

    #238U
    r238UHe = zeros(size(r238U))
    @inbounds for ri = 1:length(rSteps)
        for i=1:length(alphaRadii238U)
            # Effective radial alpha deposition from 238U
            intersectiondensity!(dInt,rEdges,relVolumes,alphaRadii238U[i],rSteps[ri])
            @turbo @. r238UHe += relVolumes[ri] * dInt * r238U[ri]
        end
    end

    #235U
    r235UHe = zeros(size(r235U))
    @inbounds for ri = 1:length(rSteps)
        for i=1:length(alphaRadii235U)
            # Effective radial alpha deposition from 235U
            intersectiondensity!(dInt, rEdges,relVolumes,alphaRadii235U[i],rSteps[ri])
            @turbo @. r235UHe += relVolumes[ri] * dInt * r235U[ri]
        end
    end

    #232Th
    r232ThHe = zeros(size(r232Th))
    @inbounds for ri = 1:length(rSteps)
        for i=1:length(alphaRadii232Th)
            # Effective radial alpha deposition from 232Th
            intersectiondensity!(dInt, rEdges,relVolumes,alphaRadii232Th[i],rSteps[ri])
            @turbo @. r232ThHe += relVolumes[ri] * dInt * r232Th[ri]
        end
    end

    # Alpha decay recoil damage
    r238Udam = 8*r238U # No smoothing of alpha damage, 8 alphas per 238U
    r235Udam = 7*r235U # No smoothing of alpha damage, 7 alphas per 235U
    r232Thdam = 6*r232Th # No smoothing of alpha damage, 6 alphas per 232 Th

    # Calculate corrected alpha deposition and recoil damage each time step for each radius
    dt_2 = dt/2
    decay = Array{T}(undef, ntSteps)
    buffer = zeros(T, ntSteps, nrSteps-2)
    # Allocate deposition and damage arrays
    alphaDeposition = zeros(T, ntSteps, nrSteps-2)
    alphaDamage = zeros(T, ntSteps, nrSteps-2)
    # 238U
    @turbo @. decay = exp(λ238U*(ageSteps + dt_2)) - exp(λ238U*(ageSteps - dt_2))
    mul!(buffer, decay, r238UHe')
    @turbo @. alphaDeposition += buffer
    mul!(buffer, decay, r238Udam')
    @turbo @. alphaDamage += buffer
    # 235U
    @turbo @. decay = exp(λ235U*(ageSteps + dt_2)) - exp(λ235U*(ageSteps - dt_2))
    mul!(buffer, decay, r235UHe')
    @turbo @. alphaDeposition += buffer
    mul!(buffer, decay, r235Udam')
    @turbo @. alphaDamage += buffer
    # 232Th
    @turbo @. decay = exp(λ232Th*(ageSteps + dt_2)) - exp(λ232Th*(ageSteps - dt_2))
    mul!(buffer, decay, r232ThHe')
    @turbo @. alphaDeposition += buffer
    mul!(buffer, decay, r232Thdam')
    @turbo @. alphaDamage += buffer

    # Allocate additional variables that will be needed for Crank-Nicholson
    annealedDamage = similar(alphaDamage)
    β = Array{T}(undef, nrSteps) # First row of annealedDamage

    # Allocate arrays for diffusivities
    Dz = Array{T}(undef, ntSteps)
    DN17 = Array{T}(undef, ntSteps)

    # Allocate output matrix for all timesteps
    u = Array{T}(undef, nrSteps, ntSteps)

    # Allocate variables for tridiagonal matrix and RHS
    dl = ones(T, nrSteps-1)    # Sub-diagonal row
    d = ones(T, nrSteps)       # Diagonal
    du = ones(T, nrSteps-1)    # Supra-diagonal row
    du2 = ones(T, nrSteps-2)   # sup-sup-diagonal row for pivoting
    ipiv = Vector{LinearAlgebra.BlasInt}(undef, nrSteps) # For pivoting

    # Tridiagonal matrix for LHS of Crank-Nicholson equation with regular grid cells
    A = Tridiagonal(dl, d, du, du2)

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = Array{T}(undef, nrSteps)

    return Zircon{T}(
        dt,
        ageSteps,
        tSteps,
        ntSteps,
        dr,
        rSteps,
        rEdges,
        relVolumes,
        nrSteps,
        r238U,
        r235U,
        r232Th,
        r238UHe,
        r235UHe,
        r232ThHe,
        alphaDeposition,
        alphaDamage,
        annealedDamage,
        u,
        β,
        Dz,
        DN17,
        A,
        ipiv,
        y,
    )
end
export Zircon

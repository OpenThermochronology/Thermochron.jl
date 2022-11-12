"""
```julia
intersectionfraction(r₁, r₂, d)
```
Calculate the fraction of the surface area of a sphere s2 with radius `r₂`
that intersects the interior of a sphere s1 of radius `r₁` if the two are
separated by distance `d`. Assumptions: `r₁`>0, `r₂`>0, `d`>=0
"""
function intersectionfraction(r₁::T, r₂::T, d::T) where T <: Number
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
function intersectiondensity(rEdges::Vector{T}, rVolumes::Vector{T}, ralpha::T, d::T) where T <: Number
    n = length(rEdges) - 1
    dInt = Array{Float64,1}(undef, n)
    fIntLast = intersectionfraction(rEdges[1],ralpha,d)
    for i=1:n
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
    A::Tridiagonal{T, Vector{T}}
    y::Vector{T}
end

function Zircon(r::T, dr::Number, Uppm::T, Thppm::T, dt::Number, ageSteps::AbstractVector{T}) where T<:Number

    # Temporal discretization
    tSteps = reverse(ageSteps)
    ntSteps = length(tSteps) # Number of time steps

    # Crystal size and spatial discretization
    rSteps = Array{T}(0+dr/2 : dr: r-dr/2)
    rEdges = Array{T}(0 : dr : r) # Edges of each radius element
    nrSteps = length(rSteps)+2 # number of radial grid points
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

    #238U
    r238UHe = zeros(size(r238U))
    for i=1:length(alphaRadii238U)
        for ri = 1:length(rSteps)
            # Effective radial alpha deposition from 238U
            r238UHe += relVolumes[ri] * intersectiondensity(rEdges,relVolumes,alphaRadii238U[i],rSteps[ri]) * r238U[ri]
        end
    end

    #235U
    r235UHe = zeros(size(r235U))
    for i=1:length(alphaRadii235U)
        for ri = 1:length(rSteps)
            # Effective radial alpha deposition from 235U
            r235UHe += relVolumes[ri] * intersectiondensity(rEdges,relVolumes,alphaRadii235U[i],rSteps[ri]) * r235U[ri]
        end
    end

    #232Th
    r232ThHe = zeros(size(r232Th))
    for i=1:length(alphaRadii232Th)
        for ri = 1:length(rSteps)
            # Effective radial alpha deposition from 232Th
            r232ThHe += relVolumes[ri] * intersectiondensity(rEdges,relVolumes,alphaRadii232Th[i],rSteps[ri]) * r232Th[ri]
        end
    end

    # Calculate corrected alpha deposition each time step for each radius
    alphaDeposition = (exp.(λ238U.*(ageSteps .+ dt/2)) - exp.(λ238U.*(ageSteps .- dt/2))) * (r238UHe') +
        (exp.(λ235U.*(ageSteps .+ dt/2)) - exp.(λ235U.*(ageSteps .- dt/2))) * (r235UHe') +
        (exp.(λ232Th.*(ageSteps .+ dt/2)) - exp.(λ232Th.*(ageSteps .- dt/2))) * (r232ThHe')

    r238Udam = 8*r238U # No smoothing of alpha damage, 8 alphas per 238U
    r235Udam = 7*r235U # No smoothing of alpha damage, 7 alphas per 235U
    r232Thdam = 6*r232Th # No smoothing of alpha damage, 6 alphas per 232 Th

    # Calculate corrected alpha damage at each time step for each radius
    alphaDamage = (exp.(λ238U.*(ageSteps .+ dt/2)) - exp.(λ238U.*(ageSteps .- dt/2))) * (r238Udam') +
        (exp.(λ235U.*(ageSteps .+ dt/2)) - exp.(λ235U.*(ageSteps .- dt/2))) * (r235Udam') +
        (exp.(λ232Th.*(ageSteps .+ dt/2)) - exp.(λ232Th.*(ageSteps .- dt/2))) * (r232Thdam')

    # Allocate additional variables that will be needed for Crank-Nicholson
    annealedDamage = similar(alphaDamage)
    β = Array{T}(undef, nrSteps) # First row of annealedDamage

    # Allocate output matrix for all timesteps
    u = Array{T}(undef, ntSteps, nrSteps)

    # Allocate variables for tridiagonal matrix and RHS
    dl = ones(T, nrSteps-1)    # Sub-diagonal row
    d = ones(T, nrSteps)       # Diagonal
    du = ones(T, nrSteps-1)    # Supra-diagonal row
    du2 = ones(T, nrSteps-2)   # sup-sup-diagonal row for pivoting

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
        A,
        y,
    )
end
export Zircon

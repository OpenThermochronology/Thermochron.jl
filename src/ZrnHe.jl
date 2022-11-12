
"""
```julia
intersectionfraction(𝓇₁, 𝓇₂, d)
```
Calculate the fraction of the surface area of a sphere s2 with radius `𝓇₂`
that intersects the interior of a sphere s1 of radius `𝓇₁` if the two are
separated by distance `d`. Assumptions: `𝓇₁`>0, `𝓇₂`>0, `d`>=0
"""
function intersectionfraction(𝓇₁::T, 𝓇₂::T, d::T) where T <: Number
    dmax = 𝓇₁+𝓇₂
    dmin = abs(𝓇₁-𝓇₂)
    if d > dmax ## If separated by more than dmax, there is no intersection
        omega = zero(T)
    elseif d > dmin
        # X is the radial distance between the center of s2 and the interseciton plane
        # See e.g., http://mathworld.wolfram.com/Sphere-SphereIntersection.html
        # x = (d^2 - 𝓇₁^2 + 𝓇₂^2) / (2 * d)

        # Let omega be is the solid angle of intersection normalized by 4pi,
        # where the solid angle of a cone is 2pi*(1-cos(theta)) and cos(theta)
        # is adjacent/hypotenuse = x/𝓇₂.
        # omega = (1 - x/𝓇₂)/2

        # Rearranged for optimization:
        omega = one(T)/2 - (d^2 - 𝓇₁^2 + 𝓇₂^2) / (4 * d * 𝓇₂)

    elseif 𝓇₁<𝓇₂ # If 𝓇₁ is entirely within 𝓇₂
        omega = zero(T)
    else
        omega = one(T) # If 𝓇₂ is entirely within 𝓇₁
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


"""
```julia
ρᵣ = DamageAnnealing(dt::Number, tSteps::Vector, TSteps::Matrix)
```
Zircon damage annealing model as in Guenthner et al. 2013 (AJS)
"""
function DamageAnnealing(dt::Number, tSteps::DenseVector, TSteps::DenseVector)
    # Allocate matrix to hold reduced track lengths for all previous timesteps
    ρᵣ = zeros(length(tSteps),length(tSteps))
    # In=-place version
    DamageAnnealing!(ρᵣ, dt, tSteps, TSteps)
end
export DamageAnnealing

"""
```julia
DamageAnnealing!(ρᵣ::Matrix, dt::Number, tSteps::Vector, TSteps::Vector)
```
In-place version of `DamageAnnealing`
"""
function DamageAnnealing!(ρᵣ::DenseMatrix, dt::Number, tSteps::DenseVector, TSteps::DenseVector)
    # Annealing model constants
    B=-0.05721
    C0=6.24534
    C1=-0.11977
    C2=-314.937 - log(1E6*365.25*24*3600) # Convert from seconds to Myr
    C3=-14.2868

    ntSteps = length(tSteps)

    # First timestep
    Teq = zeros(1,ntSteps)
    ρᵣ[1,1] = 1 / ((C0 + C1*(log(dt)-C2)/(log(1 / (TSteps[1]+273.15))-C3))^(1/B)+1)

    # All subsequent timesteps
    @inbounds for i=2:ntSteps
        # Convert any existing track length reduction for damage from
        # all previous timestep to an equivalent annealing time at the
        # current temperature
        @views Teq[1:i-1] .= exp.(C2 .+ (log(1 / (TSteps[i]+273.15)) - C3) .* (((1 ./ ρᵣ[i-1,1:i-1]) .- 1).^B .- C0) ./ C1)

        # Calculate the new reduced track lengths for all previous time steps
        # Accumulating annealing strictly in terms of reduced track length
        @views ρᵣ[i,1:i] .= 1 ./ ((C0 .+ C1 .* (log.(dt .+ Teq[1:i]) .- C2) ./ (log(1 / (TSteps[i]+273.15)) - C3)).^(1/B) .+ 1) #+273.15?

    end

    # # Guenthner et al conversion
    # map!(x->1.25*(x-0.2), ρᵣ, ρᵣ)

    # # Zero-out any reduced densities below the equivalent total annealing length
    # ρᵣ[ρᵣ.<0.36] .= 0

    # # Extrapolate from bottom of data to origin
    # ρᵣ[ρᵣ.<0.4] .= 5/8*ρᵣ[ρᵣ<0.4]

    # # Rescale reduced densities based on the equivalent total annealing length
    # ρᵣ = (ρᵣ-0.36)/(1-0.36)

    # Remove any negative reduced densities
    map!(x->ifelse(x<0., 0., x), ρᵣ, ρᵣ)

    return ρᵣ
end
export DamageAnnealing!


"""
```julia
HeAge = ZrnHeAgeSpherical(dt, ageSteps::Vector, TSteps::Vector, ρᵣ::Matrix, 𝓇, d𝓇, Uppm, Thppm)
```
Calculate the precdicted U-Th/He age of a zircon that has experienced a given t-T path
(specified by `ageSteps` for time and `TSteps` for temperature, at a time resolution of `dt`)
using a Crank-Nicholson diffusion solution for a spherical grain of radius `𝓇` at spatial resolution `d𝓇`.
"""
function ZrnHeAgeSpherical(dt::Number, ageSteps::AbstractVector{T}, TSteps::AbstractVector{T}, ρᵣ::AbstractMatrix{T}, 𝓇::T, d𝓇::Number, Uppm::T, Thppm::T, diffusionparams) where T <: Number
    # Temporal discretization
    tSteps = reverse(ageSteps)
    ntSteps = length(tSteps) # Number of time steps

    # Jaffey decay constants
    l235U = log(2)/(7.0381*10^8)*10^6 # [1/Myr]
    l238U = log(2)/(4.4683*10^9)*10^6 # [1/Myr]
    l232Th = log(2)/(1.405*10^10)*10^6 # [1/Myr]

    # Alpha stopping distances for each isotope in each decay chain, from
    # Farley et al., 1996
    alphaRadii238U = (11.78, 14.09, 13.73, 14.13, 17.32, 16.69, 28.56, 16.48,)
    alphaRadii235U = (12.58, 15.04, 19.36, 18.06, 23.07, 26.87, 22.47,)
    alphaRadii232Th = (10.99, 16.67, 18.16, 17.32, 23.61, 29.19,)

    # Other constants
    # DzEa = 165.0 # kJ/mol
    # DzD0 = 193188.0 # cm^2/sec
    # DN17Ea = 71.0 # kJ/mol
    # DN17D0 = 0.0034 #6.367E-3 # cm^2/sec
    DzEa = diffusionparams.DzEa
    DzD0 = diffusionparams.DzD0
    DN17Ea = diffusionparams.DN17Ea
    DN17D0 = diffusionparams.DN17D0


    lint0 = 45920.0 # nm
    SV = 1.669 # 1/nm
    Bα = 5.48E-19 # [g/alpha] mass of amorphous material produced per alpha decay
    Phi = 3.0
    R = 0.008314472 #kJ/(K*mol)

    # Diffusivities of crystalline and amorphous endmembers
    Dz = DzD0 .* exp.(-DzEa ./ R ./ (TSteps .+ 273.15)) # cm^2/sr
    DN17 = DN17D0 .* exp.(-DN17Ea ./ R ./ (TSteps .+ 273.15)) # cm^2/s
    Dz = Dz*10000^2*(1E6*365.25*24*3600) # Convert to micron^2/Myr
    DN17 = DN17*10000^2*(1E6*365.25*24*3600) # Convert to micron^2/Myr

    # Crystal size and spatial discretization
    rSteps = Array{Float64}(0+d𝓇/2 : d𝓇: 𝓇-d𝓇/2)
    nrSteps = length(rSteps)+2 # number of radial grid points
    rEdges = Array{Float64}(0 : d𝓇 : 𝓇) # Edges of each radius element
    relVolumes = (rEdges[2:end].^3 - rEdges[1:end-1].^3)/rEdges[end]^3 # Relative volume fraction of spherical shell corresponding to each radius element


    # Observed radial HPE profiles at present day
    r238U = Uppm.*ones(size(rSteps)) # PPM
    r235U = r238U/137.818.*ones(size(rSteps)) #PPM
    r232Th = Thppm.*ones(size(rSteps)) #PPM

    # Convert to atoms per gram
    r238U = r238U * 6.022E23 / 1E6 / 238
    r235U = r235U * 6.022E23 / 1E6 / 235
    r232Th = r232Th * 6.022E23 / 1E6 / 232


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
    alphaDeposition = (exp.(l238U.*(ageSteps .+ dt/2)) - exp.(l238U.*(ageSteps .- dt/2))) * (r238UHe') +
        (exp.(l235U.*(ageSteps .+ dt/2)) - exp.(l235U.*(ageSteps .- dt/2))) * (r235UHe') +
        (exp.(l232Th.*(ageSteps .+ dt/2)) - exp.(l232Th.*(ageSteps .- dt/2))) * (r232ThHe')

    r238Udam = 8*r238U # No smoothing of alpha damage, 8 alphas per 238U
    r235Udam = 7*r235U # No smoothing of alpha damage, 7 alphas per 235U
    r232Thdam = 6*r232Th # No smoothing of alpha damage, 6 alphas per 232 Th

    # Calculate corrected alpha damage at each time step for each radius
    alphaDamage = (exp.(l238U.*(ageSteps .+ dt/2)) - exp.(l238U.*(ageSteps .- dt/2))) * (r238Udam') +
        (exp.(l235U.*(ageSteps .+ dt/2)) - exp.(l235U.*(ageSteps .- dt/2))) * (r235Udam') +
        (exp.(l232Th.*(ageSteps .+ dt/2)) - exp.(l232Th.*(ageSteps .- dt/2))) * (r232Thdam')


    # The annealed damage matrix is the summation of the ρᵣ for each
    # previous timestep multiplied by the the alpha dose at each
    # previous timestep; this is a linear combination, which can be
    # calculated efficiently for all radii by simple matrix multiplication.
    annealedDamage = ρᵣ * alphaDamage

    # Calculate initial alpha damage
    β = Array{Float64}(undef, nrSteps)
    @turbo for k = 1:(nrSteps-2)
        fₐ = 1-exp(-Bα*annealedDamage[1,k]*Phi)
        τ = (lint0/(4.2 / ((1-exp(-Bα*annealedDamage[1,k])) * SV) - 2.5))^2
        De = 1 / ((1-fₐ)^3 / (Dz[1]/τ) + fₐ^3 / DN17[1])
        β[k+1] = 2 * d𝓇^2 / (De*dt) # Common β factor
    end
    β[1] = β[2]
    β[end] = β[end-1]

    # Allocate output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the He profile and r is radius
    u = Array{Float64}(undef, ntSteps, nrSteps)
    u[1,:] .= 0 # Initial u = v = 0 everywhere

    # Declare variables for tridiagonal matrix and RHS
    dl = ones(nrSteps-1)    # Sub-diagonal row
    d = -2.0 .- β           # Diagonal
    du = ones(nrSteps-1)    # Supra-diagonal row
    du2 = ones(nrSteps-2)   # sup-sup-diagonal row for pivoting

    y = Array{Float64}(undef, nrSteps, 1)

    # Neumann inner boundary condition (u(i,1) + u(i,2) = 0)
    d[1] = 1
    du[1] = 1

    # Dirichlet outer boundary condition (u(i,end) = u(i-1,end))
    dl[nrSteps-1] = 0
    d[nrSteps] = 1

    # Fill sparse tridiagonal matrix
    A = Tridiagonal(copy(dl), d, copy(du), du2)

    @inbounds for i=2:ntSteps

        # Calculate alpha damage
        @turbo for k = 1:(nrSteps-2)
            fₐ = 1-exp(-Bα*annealedDamage[i,k]*Phi)
            τ = (lint0/(4.2 / ((1-exp(-Bα*annealedDamage[i,k])) * SV) - 2.5))^2
            De = 1 / ((1-fₐ)^3 / (Dz[i]/τ) + fₐ^3 / DN17[i])
            β[k+1] = 2 * d𝓇^2 / (De*dt) # Common β factor
        end
        β[1] = β[2]
        β[end] = β[end-1]

        # Update tridiagonal matrix
        copyto!(A.dl, dl)       # Sub-diagonal
        @. A.d = -2. - β        # Diagonal
        copyto!(A.du, du)       # Supra-diagonal

        # Neumann inner boundary condition (u(i,1) + u(i,2) = 0)
        A.d[1] = 1
        y[1] = 0

        # Dirichlet outer boundary condition (u(i,end) = u(i-1,end))
        A.d[nrSteps] = 1
        y[nrSteps] = u[i-1,nrSteps]

        # RHS of tridiagonal Crank-Nicholson equation for regular grid cells.
        # From Ketcham, 2005
        @turbo for k = 2:nrSteps-1
            𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[i-1, k], u[i-1, k-1], u[i-1, k+1]
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊ - alphaDeposition[i, k-1]*rSteps[k-1]*β[k]
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[i,:] = A\y
        F = lu!(A)
        u[i,:] .= ldiv!(F, y)
    end

    # Convert from u (coordinate-transform'd concentration) to v (real He
    # concentration)
    vFinal = u[end,:]./[-d𝓇/2; rSteps; 𝓇+d𝓇/2]
    HeAvg = mean(vFinal[2:end-1]) # Atoms/gram

    # Raw Age (i.e., as measured)
    Avg238U = mean(r238U) # Atoms/gram
    Avg235U = mean(r235U)
    Avg232Th = mean(r232Th)
    He(t) = Avg238U*8*(exp(l238U*t)-1) + Avg235U*7*(exp(l235U*t)-1) + Avg232Th*6*(exp(l232Th*t)-1)

    # Numerically solve for helium age of the grain
    HeAge = 1
    for i=1:10
        dHe = 10000*(He(HeAge+1/10000)-He(HeAge)) # Calculate derivative
        HeAge += (HeAvg-He(HeAge))/dHe # Move towards zero (He(HeAge) == HeAvg)
    end

    return HeAge
end
export ZrnHeAgeSpherical

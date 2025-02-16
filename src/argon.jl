## -- Argon age functions

# The amount of ingrown argon since time t
calc_Ar(t, K40) = K40*BR40K*(exp(λ40K*t)-1)
# First time derivative of the amount of ingrown argon since time t
calc_dArdt(t, K40) = K40*BR40K*λ40K*exp(λ40K*t)
# Use Newton's method to solve for Ar age
function newton_ar_age(Ar::T, K40; iterations::Int=16) where {T<:Number}
    Tf = float(T)
    argonage = one(Tf)
    for _ in 1:iterations
        ∂Ar∂t = calc_dArdt(argonage, K40) # Calculate derivative
        argonage += (Ar - calc_Ar(argonage, K40))/∂Ar∂t # Move towards zero (calc_Ar(argonage) == μAr)
    end
    return max(argonage, zero(Tf))
end

## --- Calculate apparent age given a particular t-T path
"""
```julia
modelage(mineral::SphericalAr, Tsteps)
modelage(mineral::PlanarAr, Tsteps)
```
Calculate the precdicted bulk K/Ar age of a mineral that has experienced a given 
t-T path (specified by `mineral.tsteps` for time and `Tsteps` for temperature, 
at a time resolution of `step(mineral.tsteps)`) using a Crank-Nicholson diffusion 
solution for a spherical (or planar slab) grain of radius (or halfwidth ) `mineral.r` 
at spatial resolution `mineral.dr`.

Spherical implementation based on the the Crank-Nicolson solution for diffusion out of a
spherical mineral crystal in Ketcham, 2005 (doi: 10.2138/rmg.2005.58.11).
"""
function modelage(mineral::SphericalAr{T}, Tsteps::AbstractVector{T}) where T <: AbstractFloat

    # Damage and annealing constants
    D0 = mineral.D0*10000^2*SEC_MYR::T      # cm^2/sec, converted to micron^2/Myr  
    Ea = mineral.Ea::T                      # kJ/mol
    R = 0.008314472                         # kJ/(K*mol)
    ΔT = mineral.offset::T + 273.15         # Conversion from C to K, plus temperature offset from the

    # Diffusivities of crystalline and amorphous endmembers
    De = mineral.De::Vector{T}
    @assert eachindex(De) == eachindex(Tsteps)
    @turbo for i ∈ eachindex(De)
        De[i] = D0 * exp(-Ea / R / (Tsteps[i] + ΔT)) # micron^2/Myr
    end

    # Get time and radius discretization
    dr = step(mineral.rsteps)
    rsteps = mineral.rsteps
    nrsteps = mineral.nrsteps
    dt = step(mineral.tsteps)
    tsteps = mineral.tsteps
    ntsteps = length(tsteps)
    @assert eachindex(tsteps) == eachindex(Tsteps) == Base.OneTo(ntsteps)
    argondeposition = mineral.argondeposition::Matrix{T}

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}
    fill!(β, 2 * dr^2 / (De[1]*dt))

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the Ar profile and r is radius
    u = mineral.u::DenseMatrix{T}
    fill!(u, zero(T)) # initial u = v = 0 everywhere

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = mineral.y

    # Tridiagonal matrix for LHS of Crank-Nicholson equation with regular grid cells
    A = mineral.A
    fill!(A.dl, 1)          # Sub-diagonal row
    @. A.d = -2 - β         # Diagonal
    fill!(A.du, 1)          # Supra-diagonal row
    F = mineral.F               # For LU factorization

    # Neumann inner boundary condition (u[i,1] + u[i,2] = 0)
    A.d[1] = 1
    A.du[1] = 1

    # Dirichlet outer boundary condition (u[i,end] = u[i-1,end])
    A.dl[nrsteps-1] = 0
    A.d[nrsteps] = 1

    @inbounds for i = 2:ntsteps

        # Update β for current temperature
        fill!(β, 2 * dr^2 / (De[i]*dt))

        # Update tridiagonal matrix
        fill!(A.dl, 1)         # Sub-diagonal
        @. A.d = -2 - β        # Diagonal
        fill!(A.du, 1)         # Supra-diagonal

        # Neumann inner boundary condition (u(i,1) + u(i,2) = 0)
        A.du[1] = 1
        A.d[1] = 1
        y[1] = 0

        # Dirichlet outer boundary condition (u(i,end) = u(i-1,end))
        A.dl[nrsteps-1] = 0
        A.d[nrsteps] = 1
        y[nrsteps] = u[nrsteps,i-1]

        # RHS of tridiagonal Crank-Nicholson equation for regular grid cells.
        # From Ketcham, 2005 https://doi.org/10.2138/rmg.2005.58.11
        @turbo for k = 2:nrsteps-1
            𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i-1], u[k-1, i-1], u[k+1, i-1]
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊ - argondeposition[i, k-1]*rsteps[k-1]*β[k]
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i] = y
    end

    # Convert from u (coordinate-transform'd conc.) to v (real Ar conc.)
    vfinal = @views u[2:end-1,end]
    vfinal ./= rsteps
    μAr = nanmean(vfinal)                   # atoms/gram daughter
    μ40K = nanmean(mineral.r40K::Vector{T}) # atoms/gram parent

    # Numerically solve for raw Ar age of the grain (i.e., as measured)
    return newton_ar_age(μAr, μ40K)
end

function model_ll(mineral::SphericalAr, Tsteps)
    age = modelage(mineral, Tsteps)
    δ = age - mineral.age
    σ² = mineral.age_sigma^2
    -0.5*(log(2*pi*σ²) + δ^2/σ²)
end
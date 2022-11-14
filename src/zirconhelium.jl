
# Jaffey decay constants
const λ235U = log(2)/(7.0381*10^8)*10^6 # [1/Myr]
const λ238U = log(2)/(4.4683*10^9)*10^6 # [1/Myr]
const λ232Th = log(2)/(1.405*10^10)*10^6 # [1/Myr]

# The amount of ingrown helium since time t
function He(t, U238, U235, Th232)
    8*U238*(exp(λ238U*t)-1) + 7*U235*(exp(λ235U*t)-1) + 6*Th232*(exp(λ232Th*t)-1)
end
# First time derivative of the amount of ingrown helium since time t
function dHe(t, U238, U235, Th232)
    8*U238*λ238U*exp(λ238U*t) + 7*U235*λ235U*exp(λ235U*t) + 6*Th232*λ232Th*exp(λ232Th*t)
end


"""
```julia
ρᵣ = anneal(dt::Number, tSteps::Vector, TSteps::Matrix, [model::Symbol])
```
Zircon damage annealing model as in Guenthner et al. 2013 (AJS)
"""
function anneal(dt::Number, tSteps::DenseVector, TSteps::DenseVector, model=:zrdaam)
    # Allocate matrix to hold reduced track lengths for all previous timesteps
    ntSteps = length(tSteps)
    ρᵣ = zeros(ntSteps,ntSteps)
    Teq = zeros(ntSteps)
    # In=-place version
    anneal!(ρᵣ, Teq, dt, tSteps, TSteps, model)
    return ρᵣ, Teq
end
export anneal

"""
```julia
anneal!(ρᵣ::Matrix, dt::Number, tSteps::Vector, TSteps::Vector, [model::Symbol])
```
In-place version of `anneal`
"""
anneal!(ρᵣ::DenseMatrix, Teq::DenseVector, dt::Number, tSteps::DenseVector, TSteps::DenseVector, model::Symbol) = anneal!(ρᵣ, Teq, dt, tSteps, TSteps, Val(model))
function anneal!(ρᵣ::DenseMatrix, Teq::DenseVector, dt::Number, tSteps::DenseVector, TSteps::DenseVector, ::Val{:zrdaam})
    # Annealing model constants
    B=-0.05721
    C0=6.24534
    C1=-0.11977
    C2=-314.937 - log(1E6*365.25*24*3600) # Convert from seconds to Myr
    C3=-14.2868

    ntSteps = length(tSteps)
    @assert size(ρᵣ) === (ntSteps, ntSteps)
    @assert size(Teq) === (ntSteps,)
    @turbo @. Teq = 0.0

    # First timestep
    ρᵣ[1,1] = 1 / ((C0 + C1*(log(dt)-C2)/(log(1 / (TSteps[1]+273.15))-C3))^(1/B)+1)

    # All subsequent timesteps
    @inbounds for i=2:ntSteps
        lᵢ = log(1 / (TSteps[i]+273.15)) - C3

        # Convert any existing track length reduction for damage from
        # all previous timestep to an equivalent annealing time at the
        # current temperature
        pᵣᵢ = view(ρᵣ, i-1, 1:i-1)
        @turbo @. Teq[1:i-1] = exp(C2 + lᵢ * ((1/pᵣᵢ - 1)^B - C0) / C1)

        # Calculate the new reduced track lengths for all previous time steps
        # Accumulating annealing strictly in terms of reduced track length
        Teqᵢ = view(Teq, 1:i)
        @turbo @. ρᵣ[i,1:i] = 1 / ((C0 + C1 * (log(dt + Teqᵢ) - C2) / lᵢ)^(1/B) + 1)
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
export anneal!


"""
```julia
HeAge = HeAgeSpherical(zircon::Zircon, TSteps::Vector, ρᵣ::Matrix, diffusionmodel)
```
Calculate the precdicted U-Th/He age of a zircon that has experienced a given t-T path
(specified by `zircon.ageSteps` for time and `TSteps` for temperature, at a time resolution of `zircon.dt`)
using a Crank-Nicholson diffusion solution for a spherical grain of radius `zircon.r` at spatial resolution `zircon.dr`.
"""
function HeAgeSpherical(zircon::Zircon{T}, TSteps::AbstractVector{T}, ρᵣ::AbstractMatrix{T}, diffusionmodel) where T <: Number

    # Diffusion constants
    # DzEa = 165.0 # kJ/mol
    # DzD0 = 193188.0 # cm^2/sec
    # DN17Ea = 71.0 # kJ/mol
    # DN17D0 = 0.0034 #6.367E-3 # cm^2/sec
    DzEa = diffusionmodel.DzEa
    DzD0 = diffusionmodel.DzD0
    DN17Ea = diffusionmodel.DN17Ea
    DN17D0 = diffusionmodel.DN17D0

    # Damage and annealing constants
    lint0 = 45920.0 # nm
    SV = 1.669 # 1/nm
    Bα = 5.48E-19 # [g/alpha] mass of amorphous material produced per alpha decay
    Phi = 3.0
    R = 0.008314472 #kJ/(K*mol)

    # Diffusivities of crystalline and amorphous endmembers
    Dz, DN17 = zircon.Dz, zircon.DN17
    @turbo for i ∈ eachindex(Dz, DN17, TSteps)
        Dzᵢ = DzD0 * exp(-DzEa / R / (TSteps[i] + 273.15)) # cm^2/s
        Dz[i] = Dzᵢ * 10000^2*(1E6*365.25*24*3600) # Convert to micron^2/Myr
        DN17ᵢ = DN17D0 * exp(-DN17Ea / R / (TSteps[i] + 273.15)) # cm^2/s
        DN17[i] = DN17ᵢ * 10000^2*(1E6*365.25*24*3600) # Convert to micron^2/Myr
    end

    # Get time and radius discretization
    dr = zircon.dr
    rSteps = zircon.rSteps
    nrSteps = zircon.nrSteps
    dt = zircon.dt
    ntSteps = zircon.ntSteps
    alphaDeposition = zircon.alphaDeposition
    alphaDamage = zircon.alphaDamage

    # The annealed damage matrix is the summation of the ρᵣ for each
    # previous timestep multiplied by the the alpha dose at each
    # previous timestep; this is a linear combination, which can be
    # calculated efficiently for all radii by simple matrix multiplication.
    annealedDamage = zircon.annealedDamage
    mul!(annealedDamage, ρᵣ, alphaDamage)

    # Calculate initial alpha damage
    β = zircon.β
    @turbo for k = 1:(nrSteps-2)
        fₐ = 1-exp(-Bα*annealedDamage[1,k]*Phi)
        τ = (lint0/(4.2 / ((1-exp(-Bα*annealedDamage[1,k])) * SV) - 2.5))^2
        De = 1 / ((1-fₐ)^3 / (Dz[1]/τ) + fₐ^3 / DN17[1])
        β[k+1] = 2 * dr^2 / (De*dt) # Common β factor
    end
    β[1] = β[2]
    β[end] = β[end-1]

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the He profile and r is radius
    u = zircon.u
    @turbo @. u[:,1] = 0 # Initial u = v = 0 everywhere

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = zircon.y

    # Tridiagonal matrix for LHS of Crank-Nicholson equation with regular grid cells
    A = zircon.A
    @turbo @. A.dl = 1         # Sub-diagonal row
    @turbo @. A.d = -2 - β     # Diagonal
    @turbo @. A.du = 1         # Supra-diagonal row
    ipiv = zircon.ipiv         # For pivoting

    # Neumann inner boundary condition (u(i,1) + u(i,2) = 0)
    A.d[1] = 1
    A.du[1] = 1

    # Dirichlet outer boundary condition (u(i,end) = u(i-1,end))
    A.dl[nrSteps-1] = 0
    A.d[nrSteps] = 1

    @inbounds for i=2:ntSteps

        # Calculate alpha damage
        @turbo for k = 1:(nrSteps-2)
            fₐ = 1-exp(-Bα*annealedDamage[i,k]*Phi)
            τ = (lint0/(4.2 / ((1-exp(-Bα*annealedDamage[i,k])) * SV) - 2.5))^2
            De = 1 / ((1-fₐ)^3 / (Dz[i]/τ) + fₐ^3 / DN17[i])
            β[k+1] = 2 * dr^2 / (De*dt) # Common β factor
        end
        β[1] = β[2]
        β[end] = β[end-1]

        # Update tridiagonal matrix
        @turbo @. A.dl = 1         # Sub-diagonal
        @turbo @. A.d = -2 - β     # Diagonal
        @turbo @. A.du = 1         # Supra-diagonal

        # Neumann inner boundary condition (u(i,1) + u(i,2) = 0)
        A.du[1] = 1
        A.d[1] = 1
        y[1] = 0

        # Dirichlet outer boundary condition (u(i,end) = u(i-1,end))
        A.dl[nrSteps-1] = 0
        A.d[nrSteps] = 1
        y[nrSteps] = u[nrSteps,i-1]

        # RHS of tridiagonal Crank-Nicholson equation for regular grid cells.
        # From Ketcham, 2005
        @turbo for k = 2:nrSteps-1
            𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i-1], u[k-1, i-1], u[k+1, i-1]
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊ - alphaDeposition[i, k-1]*rSteps[k-1]*β[k]
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i] = A\y
        F = lu!(A, ipiv)
        ldiv!(F, y)
        @turbo @. u[:,i] = y
    end

    # Convert from u (coordinate-transform'd concentration) to v (real He
    # concentration)
    vFinal = u[:,end]./[first(rSteps)-dr; rSteps; last(rSteps)+dr]
    μHe = mean(vFinal[2:end-1]) # Atoms/gram

    # Raw Age (i.e., as measured)
    μ238U = mean(zircon.r238U) # Atoms/gram
    μ235U = mean(zircon.r235U)
    μ232Th = mean(zircon.r232Th)

    # Numerically solve for helium age of the grain
    HeAge = 1
    for i=1:10
        ∂He∂t = dHe(HeAge, μ238U, μ235U, μ232Th) # Calculate derivative
        HeAge += (μHe - He(HeAge, μ238U, μ235U, μ232Th))/∂He∂t # Move towards zero (He(HeAge) == μHe)
    end

    return HeAge
end
export HeAgeSpherical


## --- Modify LinearAlgebra.lu! to reuse du2 for pivoting

using LinearAlgebra: BlasInt, checknonsingular
# Modified from LinearAlgebra stdlib to reuse `du2` and `ipiv`
function lu!(A::Tridiagonal{T,V}, ipiv::Vector{BlasInt}, pivot::Union{RowMaximum,NoPivot} = RowMaximum(); check::Bool = true) where {T,V}

    # Extract values
    n = size(A, 1)
    @assert length(ipiv) == n

    # Initialize variables
    info = 0
    dl = A.dl
    d = A.d
    du = A.du
    if dl === du
        throw(ArgumentError("off-diagonals of `A` must not alias"))
    end
    # Check if Tridiagonal matrix already has du2 for pivoting
    has_du2_defined = isdefined(A, :du2) && isa(A.du2, V) && length(A.du2) == max(0, n-2)
    if has_du2_defined
        du2 = A.du2::V
    else
        du2 = fill!(similar(d, max(0, n-2)), 0)::V
    end


    @inbounds begin
        for i = 1:n
            ipiv[i] = i
        end
        for i = 1:n-2
            # pivot or not?
            if pivot === NoPivot() || abs(d[i]) >= abs(dl[i])
                # No interchange
                if d[i] != 0
                    fact = dl[i]/d[i]
                    dl[i] = fact
                    d[i+1] -= fact*du[i]
                    du2[i] = 0
                end
            else
                # Interchange
                fact = d[i]/dl[i]
                d[i] = dl[i]
                dl[i] = fact
                tmp = du[i]
                du[i] = d[i+1]
                d[i+1] = tmp - fact*d[i+1]
                du2[i] = du[i+1]
                du[i+1] = -fact*du[i+1]
                ipiv[i] = i+1
            end
        end
        if n > 1
            i = n-1
            if pivot === NoPivot() || abs(d[i]) >= abs(dl[i])
                if d[i] != 0
                    fact = dl[i]/d[i]
                    dl[i] = fact
                    d[i+1] -= fact*du[i]
                end
            else
                fact = d[i]/dl[i]
                d[i] = dl[i]
                dl[i] = fact
                tmp = du[i]
                du[i] = d[i+1]
                d[i+1] = tmp - fact*d[i+1]
                ipiv[i] = i+1
            end
        end
        # check for a zero on the diagonal of U
        for i = 1:n
            if d[i] == 0
                info = i
                break
            end
        end
    end
    B = has_du2_defined ? A : Tridiagonal{T,V}(dl, d, du, du2)
    check && checknonsingular(info, pivot)
    return LU{T,Tridiagonal{T,V},typeof(ipiv)}(B, ipiv, convert(BlasInt, info))
end

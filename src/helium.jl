## -- Helium age functions

# The amount of ingrown helium since time t
function calc_He(t, U238, U235, Th232, Sm147=0.0)
    8*U238*(exp(λ238U*t)-1) + 7*U235*(exp(λ235U*t)-1) + 6*Th232*(exp(λ232Th*t)-1) + Sm147*(exp(λ147Sm*t)-1)
end
# First time derivative of the amount of ingrown helium since time t
function calc_dHedt(t, U238, U235, Th232, Sm147=0.0)
    8*U238*λ238U*exp(λ238U*t) + 7*U235*λ235U*exp(λ235U*t) + 6*Th232*λ232Th*exp(λ232Th*t) + Sm147*λ147Sm*exp(λ147Sm*t)
end
# Use Newton's method to solve for He age
function newton_he_age(He::T, U238, U235, Th232, Sm147=zero(T); iterations::Int=16) where {T<:Number}
    Tf = float(T)
    heliumage = one(Tf)
    for _ in 1:iterations
        ∂He∂t = calc_dHedt(heliumage, U238, U235, Th232, Sm147) # Calculate derivative
        heliumage += (He - calc_He(heliumage, U238, U235, Th232, Sm147))/∂He∂t # Move towards zero (He(heliumage) == μHe)
    end
    return max(heliumage, zero(Tf))
end

## --- Damage and annealing functions

"""
```julia
ρᵣ = anneal(dt::Number, tsteps::Vector, Tsteps::Matrix, [model::DiffusivityModel=ZRDAAM()])
```
ZirconHe damage annealing model as in Guenthner et al. 2013 (AJS)
"""
function anneal(dt::Number, tsteps::AbstractVector, Tsteps::AbstractVector, dm::DiffusivityModel=ZRDAAM())
    # Allocate matrix to hold reduced track lengths for all previous timesteps
    ntsteps = length(tsteps)
    ρᵣ = zeros(ntsteps,ntsteps)
    teq = zeros(ntsteps)
    # In=-place version
    anneal!(ρᵣ, teq, dt, tsteps, Tsteps, dm)
    return ρᵣ, teq
end

"""
```julia
anneal!(data::Vector{<:Chronometer}, ::Type{<:HeliumSample}, tsteps, Tsteps, dm::DiffusivityModel)
anneal!(mineral::ZirconHe, Tsteps, dm::ZRDAAM)
anneal!(mineral::ApatiteHe, Tsteps, dm::RDAAM)
anneal!(ρᵣ::Matrix, dt::Number, tsteps, Tsteps, [dm::DiffusivityModel=ZRDAAM()])
```
In-place version of `anneal`
"""
function anneal!(data::Vector{<:ChronometerUnion{T}}, ::Type{C}, tsteps::AbstractRange{T}, Tsteps::AbstractVector{T}, dm::DiffusivityModel{T}) where {T<:AbstractFloat, C<:HeliumSample}
    dt = step(tsteps)::T
    if any(x->isa(x, C), data)
        imax = argmax(i->isa(data[i], C) ? length(data[i].tsteps) : 0, eachindex(data))
        m = data[imax]::C
        tmax = last(m.tsteps)::T
        @assert dt == step(m.tsteps)
        first_index = 1 + Int((last(tsteps) - tmax)÷dt)
        if first_index > 1
            anneal!(m, @views(Tsteps[first_index:end]), dm)
        else
            anneal!(m, Tsteps, dm)
        end
        pr = m.pr
        for i in eachindex(data)
            if i!=imax && isa(data[i], C)
                c = data[i]::C
                @assert dt == step(c.tsteps)
                first_index = 1 + Int((tmax - last(c.tsteps))÷dt)
                if first_index > 1
                    mul!(c.annealeddamage, @views(pr[first_index:end, first_index:end]), c.alphadamage)
                else
                    mul!(c.annealeddamage, pr, c.alphadamage)
                end
            end
        end
    end
    return data
end
function anneal!(mineral::ZirconHe, Tsteps::AbstractVector, dm::ZRDAAM)
    anneal!(mineral.pr, view(mineral.annealeddamage,:,1), step(mineral.tsteps), mineral.tsteps, Tsteps, dm)
    mul!(mineral.annealeddamage, mineral.pr, mineral.alphadamage)
    return mineral
end
function anneal!(mineral::ApatiteHe, Tsteps::AbstractVector, dm::RDAAM)
    anneal!(mineral.pr, view(mineral.annealeddamage,:,1), step(mineral.tsteps), mineral.tsteps, Tsteps, dm)
    mul!(mineral.annealeddamage, mineral.pr, mineral.alphadamage)
    return mineral
end
function anneal!(ρᵣ::AbstractMatrix{T}, teq::AbstractVector{T}, dt::Number, tsteps::AbstractVector, Tsteps::AbstractVector, dm::ZRDAAM{T}) where T <: AbstractFloat
    @assert eachindex(tsteps) == eachindex(Tsteps) == eachindex(teq) == axes(ρᵣ, 1) == axes(ρᵣ,2) == Base.OneTo(length(tsteps))
    ntsteps = length(tsteps)
    ∅ = zero(T)
    fill!(teq, ∅)

    # First timestep
    ρᵣ[1,1] = 1 / ((dm.C0 + dm.C1*(log(dt)-dm.C2)/(log(1 / (Tsteps[1]+273.15))-dm.C3))^(1/dm.beta)+1)

    # All subsequent timesteps
    @inbounds for i=2:ntsteps
        lᵢ = log(1 / (Tsteps[i]+273.15)) - dm.C3

        # Convert any existing track length reduction for damage from
        # all previous timestep to an equivalent annealing time at the
        # current temperature
        @turbo for j in 1:i-1
            teq[j] = exp(dm.C2 + lᵢ * ((1/ρᵣ[i-1, j] - 1)^dm.beta - dm.C0) / dm.C1)
        end

        # Calculate the new reduced track lengths for all previous time steps
        # Accumulating annealing strictly in terms of reduced track length
        teqᵢ = view(teq, 1:i)
        @turbo for j in 1:i
            ρᵣ[i,j] = 1 / ((dm.C0 + dm.C1 * (log(dt + teq[j]) - dm.C2) / lᵢ)^(1/dm.beta) + 1)
        end
    end

    # Guenthner et al volume-length conversion
    rmin = dm.rmin
    scale = 1/(1-rmin)
    @fastmath @inbounds for j ∈ 1:ntsteps
        for i ∈ j:ntsteps
            ρᵣ[i,j] = if ρᵣ[i,j] >= rmin
                (ρᵣ[i,j] - rmin) * scale
            else
                zero(T)
            end
        end
    end

    # # Alternative conversion: Extrapolate from bottom of data to origin
    # map!(x->(x<0.4 ? 5/8x : x), ρᵣ, ρᵣ)

    return ρᵣ
end
function anneal!(ρᵣ::AbstractMatrix{T}, teq::AbstractVector{T}, dt::Number, tsteps::AbstractVector, Tsteps::AbstractVector, dm::RDAAM{T}) where T <: AbstractFloat

    ∅ = zero(T)
    ntsteps = length(tsteps)
    @assert size(ρᵣ) === (ntsteps, ntsteps)
    @assert size(teq) === (ntsteps,)
    fill!(teq, ∅)

    # First timestep
    ρᵣ[1,1] = 1 / ((dm.C0 + dm.C1*(log(dt)-dm.C2)/(log(1 / (Tsteps[1]+273.15))-dm.C3))^(1/dm.beta)+1)

    # All subsequent timesteps
    @inbounds for i=2:ntsteps
        lᵢ = log(1 / (Tsteps[i]+273.15)) - dm.C3

        # Convert any existing track length reduction for ρᵣ from
        # all previous timestep to an equivalent annealing time at the
        # current temperature
        @turbo for j in 1:i-1
            teq[j] = exp(dm.C2 + lᵢ * ((1/ρᵣ[i-1, j] - 1)^dm.beta - dm.C0) / dm.C1)
        end

        # Calculate the new reduced track lengths for all previous time steps
        # Accumulating annealing strictly in terms of reduced track length
        @turbo for j in 1:i
            ρᵣ[i,j] = 1 / ((dm.C0 + dm.C1 * (log(dt + teq[j]) - dm.C2) / lᵢ)^(1/dm.beta) + 1)
        end
    end

    # Corrections to ρᵣ 
    rmr0 = dm.rmr0
    kappa = dm.kappa
    scale = 1/(1-rmr0)
    @fastmath @inbounds for j ∈ 1:ntsteps
        for i ∈ j:ntsteps
            # rmr0 correction
            ρᵣ[i,j] = if ρᵣ[i,j] >= rmr0
                ((ρᵣ[i,j]-rmr0)*scale)^kappa
            else
                zero(T)
            end
        end
    end

    # Convert from length to density
    @inbounds for j ∈ 1:ntsteps
        for i ∈ j:ntsteps
            ρᵣ[i,j] = reltrackdensityap(ρᵣ[i,j])
        end
    end

    return ρᵣ
end

## --- Calculate apparent age given a particular t-T path
"""
```julia
modelage(mineral::ZirconHe, Tsteps, [ρᵣ], dm::ZRDAAM)
modelage(mineral::ApatiteHe, Tsteps, [ρᵣ], dm::RDAAM)
modelage(mineral::SphericalHe, Tsteps)
modelage(mineral::PlanarHe, Tsteps)
```
Calculate the predicted bulk U-Th/He age of a zircon, apatite, or other mineral
that has experienced a given t-T path (specified by `mineral.tsteps` for time
and `Tsteps` for temperature, at a time resolution of `step(mineral.tsteps)`) 
using a Crank-Nicolson diffusion solution for a spherical (or planar slab) grain
of radius (or halfwidth) `mineral.r` at spatial resolution `mineral.dr`.

Spherical implementation based on the the Crank-Nicolson solution for diffusion out of a
spherical mineral crystal in Ketcham, 2005 (doi: 10.2138/rmg.2005.58.11).
"""
function modelage(mineral::ZirconHe, Tsteps::AbstractVector, ρᵣ::AbstractMatrix, dm::ZirconHeliumModel)
    mul!(mineral.annealeddamage, ρᵣ, mineral.alphadamage)
    modelage(mineral, Tsteps, dm)
end
function modelage(mineral::ApatiteHe, Tsteps::AbstractVector, ρᵣ::AbstractMatrix, dm::ApatiteHeliumModel)
    mul!(mineral.annealeddamage, ρᵣ, mineral.alphadamage)
    modelage(mineral, Tsteps, dm)
end
function modelage(zircon::ZirconHe{T}, Tsteps::AbstractVector{T}, dm::ZRDAAM{T}) where T <: AbstractFloat

    # Damage and annealing constants
    DzEa = dm.DzEa::T                           # kJ/mol
    DzD0 = dm.DzD0*10000^2*SEC_MYR::T           # cm^2/sec, converted to micron^2/Myr  
    DN17Ea = dm.DN17Ea::T                       # kJ/mol
    DN17D0 = dm.DN17D0*10000^2*SEC_MYR::T       # cm^2/sec, converted to micron^2/Myr  
    lint0 = dm.lint0::T                         # nm
    SV = dm.SV::T                               # 1/nm
    Bα = dm.Bα::T                               # [g/alpha] mass of amorphous material produced per alpha decay
    Phi = dm.Phi::T                             # unitless
    R = 0.008314472                             # kJ/(K*mol)
    ΔT = zircon.offset::T + 273.15              # Conversion from C to K, plus temperature offset from the

    # Diffusivities of crystalline and amorphous endmembers
    Dz = zircon.Dz::Vector{T}
    DN17 = zircon.DN17::Vector{T}
    @assert eachindex(Dz) == eachindex(DN17) == eachindex(Tsteps)
    @turbo for i ∈ eachindex(Dz)
        Dz[i] = DzD0 * exp(-DzEa / R / (Tsteps[i] + ΔT)) # micron^2/Myr
        DN17[i] = DN17D0 * exp(-DN17Ea / R / (Tsteps[i] + ΔT)) # micron^2/Myr
    end

    # Get time and radius discretization
    dr = step(zircon.rsteps)
    rsteps = zircon.rsteps
    nrsteps = zircon.nrsteps
    dt = step(zircon.tsteps)
    tsteps = zircon.tsteps
    ntsteps = length(tsteps)
    @assert eachindex(tsteps) == eachindex(Tsteps) == Base.OneTo(ntsteps)
    alphadeposition = zircon.alphadeposition::Matrix{T}

    # The annealed damage matrix is the summation of the ρᵣ for each
    # previous timestep multiplied by the the alpha dose at each
    # previous timestep; this is a linear combination, which can be
    # calculated efficiently for all radii by simple matrix multiplication.
    annealeddamage = zircon.annealeddamage::Matrix{T}

    # Calculate initial alpha damage
    β = zircon.β::Vector{T}
    @turbo for k = 1:(nrsteps-2)
        fₐ = 1-exp(-Bα*annealeddamage[1,k]*Phi)
        τ = (lint0/(4.2 / ((1-exp(-Bα*annealeddamage[1,k])) * SV) - 2.5))^2
        De = 1 / ((1-fₐ)^3 / (Dz[1]/τ) + fₐ^3 / DN17[1])
        β[k+1] = 2 * dr^2 / (De*dt) # Common β factor
    end
    β[1] = β[2]
    β[end] = β[end-1]

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the He profile and r is radius
    u = zircon.u::DenseMatrix{T}
    @assert axes(u,2) == 1:ntsteps+1
    @assert axes(u,1) == 1:nrsteps
    fill!(u, zero(T)) # initial u = v = 0 everywhere

    # Vector for RHS of Crank-Nicolson equation with regular grid cells
    y = zircon.y

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = zircon.A        # Tridiagonal matrix
    F = zircon.F        # LU object for in-place lu factorization

    @inbounds for i=2:ntsteps

        # Calculate alpha damage
        @turbo for k = 1:(nrsteps-2)
            fₐ = 1-exp(-Bα*annealeddamage[i,k]*Phi)
            τ = (lint0/(4.2 / ((1-exp(-Bα*annealeddamage[i,k])) * SV) - 2.5))^2
            De = 1 / ((1-fₐ)^3 / (Dz[i]/τ) + fₐ^3 / DN17[i])
            β[k+1] = 2 * dr^2 / (De*dt) # Common β factor
        end
        β[1] = β[2]
        β[end] = β[end-1]

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
        y[nrsteps] = u[nrsteps,i]

        # RHS of tridiagonal Crank-Nicolson equation for regular grid cells.
        # From Ketcham, 2005 https://doi.org/10.2138/rmg.2005.58.11
        @turbo for k = 2:nrsteps-1
            𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊ - alphadeposition[i, k-1]*rsteps[k-1]*β[k]
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i+1] = y
    end

    # Convert from u (coordinate-transform'd conc.) to v (real He conc.)
    vfinal = @views u[2:end-1,end]
    vfinal ./= rsteps
    μHe = nanmean(vfinal, zircon.relvolumes) # Atoms/gram

    # Parent concentrations
    μ238U = nanmean(zircon.r238U, zircon.relvolumes) # Atoms/gram
    μ235U = nanmean(zircon.r235U, zircon.relvolumes)
    μ232Th = nanmean(zircon.r232Th, zircon.relvolumes)
    μ147Sm = nanmean(zircon.r147Sm, zircon.relvolumes)

    # Numerically solve for raw helium age of the grain (i.e, as measured)
    return newton_he_age(μHe, μ238U, μ235U, μ232Th, μ147Sm)
end
function modelage(apatite::ApatiteHe{T}, Tsteps::AbstractVector{T}, dm::RDAAM{T}) where T <: AbstractFloat

    # Damage and annealing constants
    D0L = dm.D0L*10000^2*SEC_MYR::T         # cm^2/sec, converted to micron^2/Myr  
    EaL = dm.EaL::T                         # kJ/mol
    EaTrap = dm.EaTrap::T                   # kJ/mol
    etaq = dm.etaq::T                       # Durango ηq
    psi = dm.psi::T                         # unitless
    omega = dm.omega::T                     # unitless
    rhoap = dm.rhoap::T                     # g/cm^3
    L = dm.L::T                             # cm
    lambdaf = dm.lambdaf::T                 # 1/time
    lambdaD = dm.lambdaD::T                 # 1/time
    R = 0.008314472                         # kJ/(K*mol)
    ΔT = apatite.offset::T + 273.15         # Conversion from C to K, plus temperature offset from the

    # Conversion factor from alphas/g to track length cm/cm^3
    damage_conversion = rhoap*(lambdaf/lambdaD)*etaq*L

    # Diffusivities of crystalline and amorphous endmembers
    DL = apatite.DL::Vector{T}
    Dtrap = apatite.Dtrap::Vector{T}
    @assert eachindex(DL) == eachindex(Dtrap) == eachindex(Tsteps)
    @turbo for i ∈ eachindex(DL)
        DL[i] = D0L * exp(-EaL / R / (Tsteps[i] + ΔT)) # micron^2/Myr
        Dtrap[i] = exp(-EaTrap / R / (Tsteps[i] + ΔT)) # unitless
    end

    # Get time and radius discretization
    dr = step(apatite.rsteps)
    rsteps = apatite.rsteps
    nrsteps = apatite.nrsteps
    dt = step(apatite.tsteps)
    tsteps = apatite.tsteps
    ntsteps = length(tsteps)
    @assert eachindex(tsteps) == eachindex(Tsteps) == Base.OneTo(ntsteps)
    alphadeposition = apatite.alphadeposition::Matrix{T}

    # The annealed damage matrix is the summation of the ρᵣ for each
    # previous timestep multiplied by the the alpha dose at each
    # previous timestep; this is a linear combination, which can be
    # calculated efficiently for all radii by simple matrix multiplication.
    annealeddamage = apatite.annealeddamage::Matrix{T}

    # Calculate initial alpha damage
    β = apatite.β::Vector{T}
    @turbo for k = 1:(nrsteps-2)
        track_density = annealeddamage[1,k]*damage_conversion # cm/cm3
        trapDiff = psi*track_density + omega*track_density^3
        De = DL[1]/(trapDiff*Dtrap[1]+1) # micron^2/Myr
        β[k+1] = 2 * dr^2 / (De*dt) # Common β factor
    end
    β[1] = β[2]
    β[end] = β[end-1]

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the He profile and r is radius
    u = apatite.u::DenseMatrix{T}
    @assert axes(u,2) == 1:ntsteps+1
    @assert axes(u,1) == 1:nrsteps
    fill!(u, zero(T)) # initial u = v = 0 everywhere

    # Vector for RHS of Crank-Nicolson equation with regular grid cells
    y = apatite.y

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = apatite.A       # Tridiagonal matrix
    F = apatite.F       # LU object for in-place lu factorization

    @inbounds for i in Base.OneTo(ntsteps)

        # Calculate alpha damage
        @turbo for k = 1:(nrsteps-2)
            track_density = annealeddamage[i,k]*damage_conversion # cm/cm3
            trapDiff = psi*track_density + omega*track_density^3
            De = DL[i]/(trapDiff*Dtrap[i]+1) # micron^2/Myr
            β[k+1] = 2 * dr^2 / (De*dt) # Common β factor
        end
        β[1] = β[2]
        β[end] = β[end-1]

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
        y[nrsteps] = u[nrsteps,i]

        # RHS of tridiagonal Crank-Nicolson equation for regular grid cells.
        # From Ketcham, 2005 https://doi.org/10.2138/rmg.2005.58.11
        @turbo for k = 2:nrsteps-1
            𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊ - alphadeposition[i, k-1]*rsteps[k-1]*β[k]
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i+1] = y
    end

    # Convert from u (coordinate-transform'd conc.) to v (real He conc.)
    vfinal = @views u[2:end-1,end]
    vfinal ./= rsteps
    μHe = nanmean(vfinal, apatite.relvolumes) # Atoms/gram daughter

    # Parent concentrations in atoms/gram
    μ238U = nanmean(apatite.r238U, apatite.relvolumes)
    μ235U = nanmean(apatite.r235U, apatite.relvolumes)
    μ232Th = nanmean(apatite.r232Th, apatite.relvolumes)
    μ147Sm = nanmean(apatite.r147Sm, apatite.relvolumes)

    # Numerically solve for raw helium age of the grain (i.e, as measured)
    return newton_he_age(μHe, μ238U, μ235U, μ232Th, μ147Sm)
end
function modelage(mineral::SphericalHe{T}, Tsteps::AbstractVector{T}) where T <: AbstractFloat

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
    ntsteps = length(mineral.tsteps)
    alphadeposition = mineral.alphadeposition::Matrix{T}

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}
    fill!(β, 2 * dr^2 / (De[1]*dt))

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the He profile and r is radius
    u = mineral.u::DenseMatrix{T}
    @assert axes(u,2) == 1:ntsteps+1
    @assert axes(u,1) == 1:nrsteps
    fill!(u, zero(T)) # initial u = v = 0 everywhere

    # Vector for RHS of Crank-Nicolson equation with regular grid cells
    y = mineral.y

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = mineral.A       # Tridiagonal matrix
    F = mineral.F       # LU object for in-place lu factorization

    @inbounds for i in Base.OneTo(ntsteps)

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
        y[nrsteps] = u[nrsteps,i]

        # RHS of tridiagonal Crank-Nicolson equation for regular grid cells.
        # From Ketcham, 2005 https://doi.org/10.2138/rmg.2005.58.11
        @turbo for k = 2:nrsteps-1
            𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊ - alphadeposition[i, k-1]*rsteps[k-1]*β[k]
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i+1] = y
    end

    # Convert from u (coordinate-transform'd conc.) to v (real He conc.)
    vfinal = @views u[2:end-1,end]
    vfinal ./= rsteps
    μHe = nanmean(vfinal, mineral.relvolumes) # Atoms/gram

    # Parent concentrations
    μ238U = nanmean(mineral.r238U, mineral.relvolumes) # Atoms/gram
    μ235U = nanmean(mineral.r235U, mineral.relvolumes)
    μ232Th = nanmean(mineral.r232Th, mineral.relvolumes)
    μ147Sm = nanmean(mineral.r147Sm, mineral.relvolumes)

    # Numerically solve for raw helium age of the grain (i.e, as measured)
    return newton_he_age(μHe, μ238U, μ235U, μ232Th, μ147Sm)
end
function modelage(mineral::PlanarHe{T}, Tsteps::AbstractVector{T}) where T <: AbstractFloat

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
    nrsteps = mineral.nrsteps
    dt = step(mineral.tsteps)
    ntsteps = length(mineral.tsteps)
    alphadeposition = mineral.alphadeposition::Matrix{T}

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}
    fill!(β, 2 * dr^2 / (De[1]*dt))

    # Output matrix for all timesteps
    u = mineral.u::DenseMatrix{T}
    @assert axes(u,2) == 1:ntsteps+1
    @assert axes(u,1) == 1:nrsteps
    fill!(u, zero(T)) # initial u = 0 everywhere

    # Vector for RHS of Crank-Nicolson equation with regular grid cells
    y = mineral.y

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = mineral.A       # Tridiagonal matrix
    F = mineral.F       # LU object for in-place lu factorization

    @inbounds for i in Base.OneTo(ntsteps)

        # Update β for current temperature
        fill!(β, 2 * dr^2 / (De[i]*dt))

        # Update tridiagonal matrix
        fill!(A.dl, 1)         # Sub-diagonal
        @. A.d = -2 - β        # Diagonal
        fill!(A.du, 1)         # Supra-diagonal

        # Neumann inner boundary condition (-u(i,1) + u(i,2) = 0)
        A.du[1] = -1
        A.d[1] = 1
        y[1] = 0

        # Dirichlet outer boundary condition (u(i,end) = u(i-1,end))
        A.dl[nrsteps-1] = 0
        A.d[nrsteps] = 1
        y[nrsteps] = u[nrsteps,i]

        # RHS of tridiagonal Crank-Nicolson equation for regular grid cells.
        # From Ketcham, 2005 https://doi.org/10.2138/rmg.2005.58.11
        @turbo for k = 2:nrsteps-1
            𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊ - alphadeposition[i, k-1]*β[k]
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i+1] = y
    end

    # Resulting He concentration
    vfinal = @views u[2:end-1,end]
    μHe = nanmean(vfinal) # Atoms/gram

    # Parent concentrations
    μ238U = nanmean(mineral.r238U::Vector{T}) # Atoms/gram
    μ235U = nanmean(mineral.r235U::Vector{T})
    μ232Th = nanmean(mineral.r232Th::Vector{T})
    μ147Sm = nanmean(mineral.r147Sm::Vector{T})

    # Numerically solve for raw helium age of the grain (i.e, as measured)
    return newton_he_age(μHe, μ238U, μ235U, μ232Th, μ147Sm)
end

# Log likelihood for model ages
function model_ll(mineral::Union{ZirconHe,ApatiteHe}, Tsteps, dm::DiffusivityModel)
    anneal!(mineral, Tsteps, dm)
    age = modelage(mineral, Tsteps, dm)
    δ = age - mineral.age
    σ² = mineral.age_sigma^2
    -0.5*(log(2*pi*σ²) + δ^2/σ²)
end
function model_ll(mineral::Union{SphericalHe,PlanarHe}, Tsteps)
    age = modelage(mineral, Tsteps)
    δ = age - mineral.age
    σ² = mineral.age_sigma^2
    -0.5*(log(2*pi*σ²) + δ^2/σ²)
end
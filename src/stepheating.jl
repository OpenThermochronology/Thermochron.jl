## -- Functions for modelling experimental degassing schedules

function degas_initialized!(mineral::Union{PlanarHe{T}, PlanarAr{T}}, step_diffusant::AbstractVector{T}, tsteps_degassing::AbstractVector{T}, Tsteps_degassing::AbstractVector{T}, dm::Diffusivity{T}; fuse::Bool=true, diffusivityratio=one(T)) where {T}
    # Constants
    D0 = (dm.D0*10000^2)::T                 # [micron^2/sec], converted from [cm^2/sec]
    Ea = dm.Ea::T                           # [kJ/mol]
    R = 0.008314472                         # [kJ/(K*mol)]
    ΔT = T(273.15)                          # Conversion from C to K

    # Calculate effective diffusivity at each time step
    De = mineral.De::Vector{T}
    @assert firstindex(De) == firstindex(Tsteps_degassing)
    @assert lastindex(De) >= lastindex(Tsteps_degassing)
    @turbo for i ∈ eachindex(Tsteps_degassing)
        De[i] = D0 * exp(-Ea / R / (Tsteps_degassing[i] + ΔT)) * diffusivityratio # [micron^2/sec]
    end

    # Get time and radius discretization
    rsteps = mineral.rsteps
    dr = step(rsteps)
    nrsteps = mineral.nrsteps::Int
    ntsteps = length(tsteps_degassing)
    @assert eachindex(step_diffusant) == eachindex(tsteps_degassing) == eachindex(Tsteps_degassing) == Base.OneTo(ntsteps)

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the Ar profile and r is radius
    u = mineral.u::DenseMatrix{T}
    @assert axes(u,2) == 1:ntsteps+1
    @assert axes(u,1) == 1:nrsteps

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = mineral.y

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = mineral.A       # Tridiagonal matrix
    F = mineral.F       # LU object for in-place lu factorization

    total_diffusant = diffusantᵢ₋ = nanmean(@views(u[2:end-1, 1]))
    @inbounds for i in Base.OneTo(ntsteps-fuse)
        # Duration of current timestep
        dt = step_at(tsteps_degassing, i)

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

        # RHS of tridiagonal Crank-Nicholson equation for regular grid cells.
        # From Ketcham, 2005 https://doi.org/10.2138/rmg.2005.58.11
        @turbo for k = 2:nrsteps-1
            𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No ingrowth on lab timescales
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i+1] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i+1] = y

        diffusantᵢ = nanmean(@views(u[2:end-1, i+1]))
        step_diffusant[i] = max(diffusantᵢ₋ - diffusantᵢ, zero(T))
        diffusantᵢ₋ = diffusantᵢ
    end
    # Degas all remaining diffusant in last step if fuse==true
    fuse && (step_diffusant[ntsteps] = max(diffusantᵢ₋, zero(T)))

    return total_diffusant
end
function degas_initialized!(mineral::Union{SphericalHe{T}, SphericalAr{T}}, step_diffusant::AbstractVector{T}, tsteps_degassing::AbstractVector{T}, Tsteps_degassing::AbstractVector{T}, dm::Diffusivity{T}; fuse::Bool=true, diffusivityratio=one(T)) where {T}
    # Constants
    D0 = (dm.D0*10000^2)::T                 # [micron^2/sec], converted from [cm^2/sec]
    Ea = dm.Ea::T                           # [kJ/mol]
    R = 0.008314472                         # [kJ/(K*mol)]
    ΔT = T(273.15)                          # Conversion from C to K

    # Calculate effective diffusivity at each time step
    De = mineral.De::Vector{T}
    @assert firstindex(De) == firstindex(Tsteps_degassing)
    @assert lastindex(De) >= lastindex(Tsteps_degassing)
    @turbo for i ∈ eachindex(Tsteps_degassing)
        De[i] = D0 * exp(-Ea / R / (Tsteps_degassing[i] + ΔT)) * diffusivityratio # [micron^2/sec]
    end

    # Get time and radius discretization
    rsteps = mineral.rsteps
    dr = step(rsteps)
    nrsteps = mineral.nrsteps::Int
    relvolumes = mineral.relvolumes::Vector{T}
    ntsteps = length(tsteps_degassing)
    @assert eachindex(step_diffusant) == eachindex(tsteps_degassing) == eachindex(Tsteps_degassing) == Base.OneTo(ntsteps)

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the Ar profile and r is radius
    u = mineral.u::DenseMatrix{T}
    @assert axes(u,2) == 1:ntsteps+1
    @assert axes(u,1) == 1:nrsteps

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = mineral.y

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = mineral.A       # Tridiagonal matrix
    F = mineral.F       # LU object for in-place lu factorization

    @inbounds for i in Base.OneTo(ntsteps-fuse)
        # Duration of current timestep
        dt = step_at(tsteps_degassing, i)

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

        # RHS of tridiagonal Crank-Nicholson equation for regular grid cells.
        # From Ketcham, 2005 https://doi.org/10.2138/rmg.2005.58.11
        @turbo for k = 2:nrsteps-1
            𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No diffusant ingrowth on lab timescales
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i+1] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i+1] = y
    end

    # Calculate diffusant lost to diffusion at each step
    total_diffusant = diffusantᵢ₋ = nanmean(@views(u[2:end-1, 1])./=rsteps, relvolumes)
    @inbounds for i in Base.OneTo(ntsteps-fuse)
        diffusantᵢ = nanmean(@views(u[2:end-1, i+1])./=rsteps, relvolumes)
        step_diffusant[i] = max(diffusantᵢ₋ - diffusantᵢ, zero(T))
        diffusantᵢ₋ = diffusantᵢ
    end
    # Degas all remaining diffusant in last step if fuse==true
    fuse && (step_diffusant[ntsteps] = max(diffusantᵢ₋, zero(T)))

    return total_diffusant
end
function degas_initialized!(mineral::ZirconHe{T}, step_diffusant::AbstractVector{T}, tsteps_degassing::AbstractVector{T}, Tsteps_degassing::AbstractVector{T}, dm::ZRDAAM{T}; fuse::Bool=true, diffusivityratio=one(T)) where {T}
    # Damage and annealing constants
    DzEa = dm.DzEa::T                           # [kJ/mol]
    DzD0 = (dm.DzD0*10000^2)::T                 # [micron^2/sec], converted from [cm^2/sec]
    DN17Ea = dm.DN17Ea::T                       # [kJ/mol]
    DN17D0 = (dm.DN17D0*10000^2)::T             # [micron^2/sec], converted from [cm^2/sec]
    lint0 = dm.lint0::T                         # [nm]
    SV = dm.SV::T                               # [1/nm]
    Bα = dm.Bα::T                               # [g/alpha] mass of amorphous material produced per alpha decay
    Phi = dm.Phi::T                             # [unitless]
    R = 0.008314472                             # [kJ/(K*mol)]
    ΔT = T(273.15)                              # Conversion from C to K

    # Diffusivities of crystalline and amorphous endmembers at each timestep
    Dz = mineral.Dz::Vector{T}
    DN17 = mineral.DN17::Vector{T}
    @assert eachindex(Dz) == eachindex(DN17) == eachindex(Tsteps_degassing)
    @turbo for i ∈ eachindex(Dz)
        Dz[i] = DzD0 * exp(-DzEa / R / (Tsteps_degassing[i] + ΔT)) # [micron^2/sec]
        DN17[i] = DN17D0 * exp(-DN17Ea / R / (Tsteps_degassing[i] + ΔT)) # [micron^2/sec]
    end

    # Get time and radius discretization
    rsteps = mineral.rsteps
    dr = step(rsteps)
    nrsteps = mineral.nrsteps::Int
    relvolumes = mineral.relvolumes::Vector{T}
    ntsteps = length(tsteps_degassing)
    @assert eachindex(step_diffusant) == eachindex(tsteps_degassing) == eachindex(Tsteps_degassing) == Base.OneTo(ntsteps)

    # The annealed damage matrix
    # We will use the last timestep of the geological annealed damage matrix,
    # Assuming He diffuses out faster than annealing occurs during laboratory degassing
    annealeddamage = mineral.annealeddamage::Matrix{T}

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the Ar profile and r is radius
    u = mineral.u::DenseMatrix{T}
    @assert axes(u,2) == 1:ntsteps+1
    @assert axes(u,1) == 1:nrsteps

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = mineral.y

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = mineral.A       # Tridiagonal matrix
    F = mineral.F       # LU object for in-place lu factorization

    @inbounds for i in Base.OneTo(ntsteps-fuse)
        # Duration of current timestep
        dt = step_at(tsteps_degassing, i)

        # Calculate alpha damage and β factor at each radius at current temperature
        @turbo for k = 1:(nrsteps-2)
            fₐ = 1-exp(-Bα*annealeddamage[ntsteps,k]*Phi)
            τ = (lint0/(4.2 / ((1-exp(-Bα*annealeddamage[ntsteps,k])) * SV) - 2.5))^2
            De = 1 / ((1-fₐ)^3 / (Dz[i]/τ) + fₐ^3 / DN17[i]) * diffusivityratio # [micron^2/sec]
            β[k+1] = 2 * dr^2 / (De*dt) # Shifted by 1 because β[1] is implicit point at negative radius
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

        # RHS of tridiagonal Crank-Nicholson equation for regular grid cells.
        # From Ketcham, 2005 https://doi.org/10.2138/rmg.2005.58.11
        @turbo for k = 2:nrsteps-1
            𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No diffusant ingrowth on lab timescales
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i+1] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i+1] = y
    end

    # Calculate diffusant lost to diffusion at each step
    total_diffusant = diffusantᵢ₋ = nanmean(@views(u[2:end-1, 1])./=rsteps, relvolumes)
    @inbounds for i in Base.OneTo(ntsteps-fuse)
        diffusantᵢ = nanmean(@views(u[2:end-1, i+1])./=rsteps, relvolumes)
        step_diffusant[i] = max(diffusantᵢ₋ - diffusantᵢ, zero(T))
        diffusantᵢ₋ = diffusantᵢ
    end
    # Degas all remaining diffusant in last step if fuse==true
    fuse && (step_diffusant[ntsteps] = max(diffusantᵢ₋, zero(T)))

    return total_diffusant
end
function degas_initialized!(mineral::ApatiteHe{T}, step_diffusant::AbstractVector{T}, tsteps_degassing::AbstractVector{T}, Tsteps_degassing::AbstractVector{T}, dm::RDAAM{T}; fuse::Bool=true, diffusivityratio=one(T)) where {T}
    # Damage and annealing constants
    D0L = (dm.D0L*10000^2)::T               # [micron^2/sec], converted from [cm^2/sec]  
    EaL = dm.EaL::T                         # [kJ/mol]
    EaTrap = dm.EaTrap::T                   # [kJ/mol]
    etaq = dm.etaq::T                       # Durango ηq
    psi = dm.psi::T                         # [unitless]
    omega = dm.omega::T                     # [unitless]
    rhoap = dm.rhoap::T                     # [g/cm^3]
    L = dm.L::T                             # [cm]
    lambdaf = dm.lambdaf::T                 # [1/time]
    lambdaD = dm.lambdaD::T                 # [1/time]
    R = 0.008314472                         # [kJ/(K*mol)]
    ΔT = T(273.15)                          # Conversion from C to K

    # Conversion factor from alphas/g to track length cm/cm^3
    damage_conversion = rhoap*(lambdaf/lambdaD)*etaq*L

    # Normal and trapping diffusivities at each timestep
    DL = mineral.DL::Vector{T}
    Dtrap = mineral.Dtrap::Vector{T}
    @assert eachindex(DL) == eachindex(Dtrap) == eachindex(Tsteps_degassing)
    @turbo for i ∈ eachindex(DL)
        DL[i] = D0L * exp(-EaL / R / (Tsteps_degassing[i] + ΔT)) # [micron^2/Myr]
        Dtrap[i] = exp(-EaTrap / R / (Tsteps_degassing[i] + ΔT)) # [unitless]
    end

    # Get time and radius discretization
    rsteps = mineral.rsteps
    dr = step(rsteps)
    nrsteps = mineral.nrsteps::Int
    relvolumes = mineral.relvolumes::Vector{T}
    ntsteps = length(tsteps_degassing)
    @assert eachindex(step_diffusant) == eachindex(tsteps_degassing) == eachindex(Tsteps_degassing) == Base.OneTo(ntsteps)

    # The annealed damage matrix
    # We will use the last timestep of the geological annealed damage matrix,
    # Assuming He diffuses out faster than annealing occurs during laboratory degassing
    annealeddamage = mineral.annealeddamage::Matrix{T}

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the Ar profile and r is radius
    u = mineral.u::DenseMatrix{T}
    @assert axes(u,2) == 1:ntsteps+1
    @assert axes(u,1) == 1:nrsteps

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = mineral.y

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = mineral.A       # Tridiagonal matrix
    F = mineral.F       # LU object for in-place lu factorization

    @inbounds for i in Base.OneTo(ntsteps-fuse)
        # Duration of current timestep
        dt = step_at(tsteps_degassing, i)

        # Calculate alpha damage and β factor at each radius at current temperature
        @turbo for k = 1:(nrsteps-2)
            track_density = annealeddamage[ntsteps,k]*damage_conversion # [cm/cm3]
            trapDiff = psi*track_density + omega*track_density^3
            De = DL[i]/(trapDiff*Dtrap[i]+1) * diffusivityratio # [micron^2/sec]
            β[k+1] = 2 * dr^2 / (De*dt) # Shifted by 1 because β[1] is implicit point at negative radius
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

        # RHS of tridiagonal Crank-Nicholson equation for regular grid cells.
        # From Ketcham, 2005 https://doi.org/10.2138/rmg.2005.58.11
        @turbo for k = 2:nrsteps-1
            𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No diffusant ingrowth on lab timescales
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i+1] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i+1] = y
    end

    # Calculate diffusant lost to diffusion at each step
    total_diffusant = diffusantᵢ₋ = nanmean(@views(u[2:end-1, 1])./=rsteps, relvolumes)
    @inbounds for i in Base.OneTo(ntsteps-fuse)
        diffusantᵢ = nanmean(@views(u[2:end-1, i+1])./=rsteps, relvolumes)
        step_diffusant[i] = max(diffusantᵢ₋ - diffusantᵢ, zero(T))
        diffusantᵢ₋ = diffusantᵢ
    end
    # Degas all remaining diffusant in last step if fuse==true
    fuse && (step_diffusant[ntsteps] = max(diffusantᵢ₋, zero(T)))

    return total_diffusant
end

# Initialize and degas daughter isotopes
function degas_daughter!(mineral::Union{HeliumSample{T}, ArgonSample{T}}, tsteps_degassing, Tsteps_degassing, dm; fuse::Bool=true) where {T}
    # Erase previous diffusion profiles
    u = fill!(mineral.u, zero(T))
    u[:,1] .= mineral.y # Initialize with final profile from prevous (geologic) inversion
    # u[end,1] = zero(T) # Zero concentration at outer boundary, since we degas in a vacuum
    # Degas and return total daughter
    return degas_initialized!(mineral, mineral.step_daughter, tsteps_degassing, Tsteps_degassing, dm; fuse)
end

# Initialize and degas tracer isotopes
function degas_tracer!(mineral::Union{PlanarAr{T}, PlanarHe{T}}, initial_tracer, tsteps_degassing, Tsteps_degassing, dm; fuse::Bool=true)  where {T}
    # Erase previous diffusion profiles
    u = fill!(mineral.u, zero(T))
    u[2:end-1,1] .= initial_tracer
    u[1,1] = u[2,1] # Symmetric inner boundary condition
    u[end,1] = zero(T) # Zero concentration at outer boundary, since we degas in a vacuum
    # Degas and return total tracer
    diffusivityratio = tracerdiffusivityratio(mineral)
    return degas_initialized!(mineral, mineral.step_tracer, tsteps_degassing, Tsteps_degassing, dm; fuse, diffusivityratio)
end
function degas_tracer!(mineral::Union{SphericalAr{T}, SphericalHe{T}, ApatiteHe{T}, ZirconHe{T}}, initial_tracer, tsteps_degassing, Tsteps_degassing, dm; fuse::Bool=true)  where {T}
    # Erase previous diffusion profiles
    u = fill!(mineral.u, zero(T))
    u[2:end-1,1] .= initial_tracer
    u[2:end-1,1] .*= mineral.rsteps # U-transform for Crank-Nicholson
    u[1,1] = -u[2,1] # Symmetric inner boundary given U-transform
    u[end,1] = zero(T) # Zero concentration at outer boundary, since we degas in a vacuum
    # Degas and return total tracer
    diffusivityratio = tracerdiffusivityratio(mineral)
    return degas_initialized!(mineral, mineral.step_tracer, tsteps_degassing, Tsteps_degassing, dm; fuse, diffusivityratio)
end
tracerdiffusivityratio(x::ArgonSample{T}) where {T} = T((40/39)^0.3)
tracerdiffusivityratio(x::HeliumSample{T}) where {T} = T((4/3)^0.3)

# Top-level degassing functions
function degas!(mineral::HeliumSample, tsteps_degassing, Tsteps_degassing, dm; fuse=true, redegastracer=false)
    total_daughter = degas_daughter!(mineral, tsteps_degassing, Tsteps_degassing, dm; fuse)
    # Now diffuse parent isotope tracer (He-3), if neccesary
    if redegastracer || !(0 < sum(mineral.step_tracer))
        degas_tracer!(mineral, total_daughter, tsteps_degassing, Tsteps_degassing, dm; fuse)
    end
    return mineral.step_tracer, mineral.step_daughter
end
function degas!(mineral::ArgonSample, tsteps_degassing, Tsteps_degassing, dm; fuse=true, redegastracer=false)
    degas_daughter!(mineral, tsteps_degassing, Tsteps_degassing, dm; fuse)
    # Now diffuse parent isotope tracer (Ar-39), if neccesary
    if redegastracer || !(0 < sum(mineral.step_tracer))
        degas_tracer!(mineral, mineral.r40K, tsteps_degassing, Tsteps_degassing, dm; fuse)
    end
    return mineral.step_tracer, mineral.step_daughter
end


## --- Age and likelihood functions for step heating data

function modelage(sdd::SingleDomain{T,<:HeliumSample}, Tsteps::AbstractVector, dm::DiffusivityModel{T}; redegastracer::Bool=false) where {T<:AbstractFloat}
    stepratio = fill!(sdd.model_age, zero(T))
    fraction = fill!(sdd.model_fraction, zero(T))
    # Degas
    age = modelage(sdd.domain, Tsteps, dm)
    tracer, daughter = degas!(sdd.domain, sdd.tsteps_degassing, sdd.Tsteps_degassing, dm; sdd.fuse, redegastracer)
    # Calculate Rstep/Rbulk for each degassing step
    @inbounds for i in eachindex(tracer, daughter)
        stepratio[i] = daughter[i]/tracer[i]
    end
    # Cumulative fraction of tracer degassed
    cumsum!(fraction, tracer)
    fraction ./= last(fraction)

    return age, stepratio, fraction
end
function modelage(sdd::SingleDomain{T,<:ArgonSample}, Tsteps::AbstractVector, dm::DiffusivityModel{T}; redegastracer::Bool=false) where {T<:AbstractFloat}
    stepage = fill!(sdd.model_age, zero(T))
    fraction = fill!(sdd.model_fraction, zero(T))
    # Degas
    age = modelage(sdd.domain, Tsteps, dm)
    tracer, daughter = degas!(sdd.domain, sdd.tsteps_degassing, sdd.Tsteps_degassing, dm; sdd.fuse, redegastracer)

    # Calculate ages for each degassing step
    for i in eachindex(tracer, daughter)
        stepage[i] = newton_ar_age(daughter[i], tracer[i])
    end
    # Cumulative fraction of tracer degassed
    cumsum!(fraction, tracer)
    fraction ./= last(fraction)

    return age, stepage, fraction
end
function modelage(mdd::MultipleDomain{T,<:ArgonSample}, Tsteps::AbstractVector, dm::MDDiffusivity{T}; redegastracer::Bool=false) where {T<:AbstractFloat}
    age = fill!(mdd.model_age, zero(T))
    tracer = fill!(mdd.model_tracer, zero(T))
    daughter = fill!(mdd.model_daughter, zero(T))
    fraction = fill!(mdd.model_fraction, zero(T))
    fuse = mdd.fuse::Bool
    # Degas
    for i in eachindex(mdd.domains, mdd.volume_fraction)
        domain = mdd.domains[i]
        modelage(domain, Tsteps, dm[i])
        p, d = degas!(domain, mdd.tsteps_degassing, mdd.Tsteps_degassing, dm[i]; fuse, redegastracer)
        @. tracer += p * mdd.volume_fraction[i]
        @. daughter += d * mdd.volume_fraction[i]
    end
    # Calculate ages for each degassing step
    for i in eachindex(tracer, daughter)
        age[i] = newton_ar_age(daughter[i], tracer[i])
    end
    # Cumulative fraction of tracer degassed
    cumsum!(fraction, tracer)
    fraction ./= last(fraction)

    return age, fraction
end

function model_ll(dd::Union{SingleDomain{T},MultipleDomain{T}}, σ::T=zero(T); rescale=false) where {T<:AbstractFloat}
    ll = zero(T)
    @inbounds for i in eachindex(dd.step_age, dd.step_age_sigma, dd.midpoint_experimental, dd.fit)
        if dd.fit[i]
            model_ageᵢ = linterp1(dd.model_fraction, dd.model_age, dd.midpoint_experimental[i])
            ll += norm_ll(dd.step_age[i], dd.step_age_sigma[i], model_ageᵢ, σ)
        end
    end
    rescale && (ll /= sqrt(count(dd.fit)))
    return ll
end

function cumulative_fraction_uncertainty(sigma, i::Int)
    i₋ = firstindex(sigma)
    i₊ = lastindex(sigma)
    σ²₋ = sum(abs2, view(sigma, i₋:i))
    σ²₊ = sum(abs2, view(sigma,(i+1):i₊))
    σ = sqrt(1/(1/σ²₋ + 1/σ²₊))
end

function cumulative_degassing_ll(dd::Union{SingleDomain{T},MultipleDomain{T}}; rescale=false) where {T<:AbstractFloat}
    ll = zero(T)
    fit_until = findlast(dd.fit)
    @inbounds for i in eachindex(dd.tsteps_experimental, dd.fraction_experimental, dd.fraction_experimental_sigma, dd.fit)
        if i <= fit_until
            σ = cumulative_fraction_uncertainty(dd.fraction_experimental_sigma, i)
            σ ≈ 0 && continue
            model_fractionᵢ = linterp1(dd.tsteps_degassing, dd.model_fraction, dd.tsteps_experimental[i])
            ll += norm_ll(dd.fraction_experimental[i], σ, model_fractionᵢ)
        end
    end
    rescale && (ll /= sqrt(fit_until))
    return ll
end

function stepwise_degassing_ll(dd::Union{SingleDomain{T},MultipleDomain{T}}; rescale=false) where {T<:AbstractFloat}
    ll = zero(T)
    last_model_fractionᵢ = zero(T)
    last_fraction_experimentalᵢ = zero(T)
    @inbounds for i in eachindex(dd.tsteps_experimental, dd.fraction_experimental, dd.fraction_experimental_sigma, dd.fit)
        model_fractionᵢ = linterp1(dd.tsteps_degassing, dd.model_fraction, dd.tsteps_experimental[i])
        if dd.fit[i]
            δmodel = model_fractionᵢ - last_model_fractionᵢ
            δexperimental = dd.fraction_experimental[i] - last_fraction_experimentalᵢ
            ll += norm_ll(δexperimental, dd.fraction_experimental_sigma[i], δmodel)
        end
        last_model_fractionᵢ₋ = model_fractionᵢ
        last_fraction_experimentalᵢ₋ = dd.fraction_experimental[i]
    end
    rescale && (ll /= sqrt(count(dd.fit)))
    return ll
end

## --- End of File
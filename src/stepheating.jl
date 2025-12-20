## -- Functions for modelling experimental degassing schedules

const PlanarNobleGas{T} = Union{PlanarHe{T}, PlanarAr{T}}
const SphericalNobleGas{T} = Union{ZirconHe{T},ApatiteHe{T},SphericalHe{T},SphericalAr{T}}

function updatebeta!(β::Vector{T}, mineral::ZirconHe{T}, dm::ZRDAAM{T}, dt::T, TK::T, tstep::Int, diffusivityratio=one(T); setting::Symbol=:geological) where {T}
    # The annealed damage matrix
    # We will use the last timestep of the geological annealed damage matrix,
    # Assuming He diffuses out faster than annealing occurs during laboratory degassing
    annealeddamage = mineral.annealeddamage::Matrix{T}

    # Damage and annealing constants
    DzEa = dm.DzEa::T                           # [kJ/mol]
    DN17Ea = dm.DN17Ea::T                       # [kJ/mol]
    lint0 = dm.lint0::T                         # [nm]
    SV = dm.SV::T                               # [1/nm]
    Bα = dm.Bα::T                               # [g/alpha] mass of amorphous material produced per alpha decay
    Phi = dm.Phi::T                             # [unitless]
    R = 0.008314472                             # [kJ/(K*mol)]
    if setting === :laboratory
        DzD0 = (dm.DzD0*10000^2)::T             # [micron^2/sec], converted from [cm^2/sec]
        DN17D0 = (dm.DN17D0*10000^2)::T         # [micron^2/sec], converted from [cm^2/sec]
        damagestep = lastindex(annealeddamage, 1)
    else
        DzD0 = (dm.DzD0*10000^2*SEC_MYR)::T     # [micron^2/Myr], converted from [cm^2/sec]
        DN17D0 = (dm.DN17D0*10000^2*SEC_MYR)::T # [micron^2/Myr], converted from [cm^2/sec]
        damagestep = tstep
    end

    # Radial step size
    dr = step(mineral.rsteps)::T

    # Endmember diffusion constants
    Dz = DzD0 * exp(-DzEa / (R * TK)) * diffusivityratio # [micron^2/sec]
    DN17 = DN17D0 * exp(-DN17Ea / (R * TK)) * diffusivityratio # [micron^2/sec

    # Each radial step except first and latst
    @inbounds @simd ivdep for k in (eachindex(β)[2:end-1])
        dam = annealeddamage[damagestep,k-1] # Shifted by 1 because β[1] is implicit point at negative radius
        fₐ = 1-exp(-Bα*dam*Phi)
        τ = (lint0/(4.2 / ((1-exp(-Bα*dam)) * SV) - 2.5))^2
        De = 1 / ((1-fₐ)^3 / (Dz/τ) + fₐ^3 / DN17)  # Effetive diffusivity
        β[k] = 2 * dr^2 / (De*dt)
    end
    # First and last radial step
    β[1] = β[2]
    β[end] = β[end-1]

    return β
end
function updatebeta!(β::Vector{T}, mineral::ApatiteHe{T}, dm::RDAAM{T}, dt::T, TK::T, tstep::Int, diffusivityratio=one(T); setting::Symbol=:laboratory) where {T}
    # The annealed damage matrix
    # We will use the last timestep of the geological annealed damage matrix,
    # Assuming He diffuses out faster than annealing occurs during laboratory degassing
    annealeddamage = mineral.annealeddamage::Matrix{T}

    # Damage and annealing constants
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
    if setting === :laboratory
        D0L = (dm.D0L*10000^2)::T           # [micron^2/sec], converted from [cm^2/sec]  
        damagestep = lastindex(annealeddamage, 1)
    else
        D0L = (dm.D0L*10000^2*SEC_MYR)::T   # [micron^2/Myr], converted from [cm^2/sec]  
        damagestep = tstep
    end

    # Radial step size
    dr = step(mineral.rsteps)::T

    # Conversion factor from alphas/g to track length cm/cm^3
    damage_conversion = rhoap*(lambdaf/lambdaD)*etaq*L

    # Normal and trapping diffusivities
    DL = D0L * exp(-EaL / (R * TK)) * diffusivityratio # [micron^2/t]
    Dtrap = exp( EaTrap / (R * TK)) # [unitless]

    # Each radial step except first and latst
    @inbounds @simd ivdep for k in (eachindex(β)[2:end-1])
        track_density = annealeddamage[damagestep, k-1]*damage_conversion # [cm/cm3]
        trap = (psi*track_density + omega*track_density^3)*Dtrap
        De = DL/(trap+1) # [micron^2/t]
        β[k] = 2 * dr^2 / (De*dt)
    end
    # First and last radial step
    β[1] = β[2]
    β[end] = β[end-1]

    return β
end
function updatebeta!(β::Vector{T}, mineral::NobleGasSample{T}, dm::Diffusivity{T}, dt::T, TK::T, ::Int, diffusivityratio=one(T); setting::Symbol=:laboratory) where {T}
    # Constants
    Ea = dm.Ea::T                           # [kJ/mol]
    R = 0.008314472                         # [kJ/(K*mol)]
    if setting === :laboratory
        D0 = (dm.D0*10000^2)::T             # [micron^2/sec], converted from [cm^2/sec]
    else
        D0 = (dm.D0*10000^2*SEC_MYR)::T     # [micron^2/Myr], converted from [cm^2/sec]
    end

    # Radial step size
    dr = step(mineral.rsteps)::T

    # Diffusivity
    De = D0 * exp(-Ea / (R * TK)) * diffusivityratio # [micron^2/t]

    # Update β for current temperature
    fill!(β, 2 * dr^2 / (De*dt))

    return β
end

function crank_nicolson!(mineral::PlanarNobleGas{T}, tsteps_degassing::AbstractVector{T}, Tsteps_degassing::AbstractVector{T}, dm::Diffusivity{T}; fuse::Bool=true, diffusivityratio=one(T), setting::Symbol=:laboratory) where {T}
    @assert setting===:laboratory || setting===:geological

    # Temperature offset
    if setting===:laboratory
        ΔT = T(273.15)                              # Conversion from C to K
    else
        ΔT = T(273.15) + mineral.offset             # Conversion from C to K
    end

    # Check time and radius discretization
    nrsteps = mineral.nrsteps::Int
    ntsteps = length(tsteps_degassing)
    @assert eachindex(tsteps_degassing) == eachindex(Tsteps_degassing) == Base.OneTo(ntsteps)

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

        # Beta factor containing diffusivity information for each radial step
        updatebeta!(β, mineral, dm, dt, Tsteps_degassing[i]+ΔT, i, diffusivityratio; setting)

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
    end

    return mineral
end
function crank_nicolson!(mineral::SphericalNobleGas{T}, tsteps_degassing::AbstractVector{T}, Tsteps_degassing::AbstractVector{T}, dm::DiffusivityModel{T}; fuse::Bool=true, diffusivityratio::Number=one(T), setting::Symbol=:laboratory) where {T}
    @assert setting===:laboratory || setting===:geological

    # Temperature offset
    if setting===:laboratory
        ΔT = T(273.15)                              # Conversion from C to K
    else
        ΔT = T(273.15) + mineral.offset             # Conversion from C to K
    end

    # Check time and radius discretization
    nrsteps = mineral.nrsteps::Int
    ntsteps = length(tsteps_degassing)
    @assert eachindex(tsteps_degassing) == eachindex(Tsteps_degassing) == Base.OneTo(ntsteps)

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

        # Beta factor containing diffusivity information for each radial step
        updatebeta!(β, mineral, dm, dt, Tsteps_degassing[i]+ΔT, i, diffusivityratio; setting)

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

    return mineral
end

function diffusant_lost!(step_diffusant::AbstractVector{T}, mineral::PlanarNobleGas{T}; fuse::Bool=true) where {T}
    # Check time and radius discretization
    ntsteps = length(step_diffusant)
    @assert eachindex(step_diffusant) == Base.OneTo(ntsteps)

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the Ar profile and r is radius
    u = mineral.u::DenseMatrix{T}
    @assert axes(u,2) == Base.OneTo(ntsteps+1)
    @assert axes(u,1) == Base.OneTo(mineral.nrsteps)

    # Calculate diffusant lost to diffusion at each step
    total_diffusant = diffusantᵢ₋ = nanmean(@views(u[2:end-1, 1]))
    @inbounds for i in Base.OneTo(ntsteps-fuse)
        diffusantᵢ = nanmean(@views(u[2:end-1, i+1]))
        step_diffusant[i] = max(diffusantᵢ₋ - diffusantᵢ, zero(T))
        diffusantᵢ₋ = diffusantᵢ
    end
    # Degas all remaining diffusant in last step if fuse==true
    fuse && (step_diffusant[ntsteps] = max(diffusantᵢ₋, zero(T)))

    return total_diffusant
end
function diffusant_lost!(step_diffusant::AbstractVector{T}, mineral::SphericalNobleGas{T}; fuse::Bool=true) where {T}
    # Check time and radius discretization
    ntsteps = length(step_diffusant)
    @assert eachindex(step_diffusant) == Base.OneTo(ntsteps)
    nrsteps = mineral.nrsteps::Int
    rsteps = mineral.rsteps::AbstractVector{T}
    relvolumes = mineral.relvolumes::Vector{T}
    @assert eachindex(relvolumes) == eachindex(rsteps) == Base.OneTo(nrsteps-2)

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the Ar profile and r is radius
    u = mineral.u::DenseMatrix{T}
    @assert axes(u,2) == Base.OneTo(ntsteps+1)
    @assert axes(u,1) == Base.OneTo(nrsteps)

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

## --- Initialize and degas daughter isotopes

function degas_daughter!(mineral::NobleGasSample{T}, tsteps_degassing, Tsteps_degassing, dm; fuse::Bool=true) where {T}
    # Erase previous diffusion profiles
    u = fill!(mineral.u, zero(T))
    u[:,1] .= mineral.y # Initialize with final profile from prevous (geologic) inversion
    # u[end,1] = zero(T) # Zero concentration at outer boundary, since we degas in a vacuum
    # Degas and return total daughter
    crank_nicolson!(mineral, tsteps_degassing, Tsteps_degassing, dm; fuse)
    return diffusant_lost!(mineral.step_daughter, mineral; fuse)
end

## --- Initialize and degas tracer isotopes
function degas_tracer!(mdd::MultipleDomain{T}, dm::MDDiffusivity{T}) where {T<:AbstractFloat}
    fraction = fill!(mdd.model_fraction, zero(T))
    tracer = fill!(mdd.model_tracer, zero(T))
    # Degas
    for i in eachindex(mdd.domains, mdd.volume_fraction)
        domain = mdd.domains[i]
        degas_tracer!(domain, one(T), mdd.tsteps_degassing, mdd.Tsteps_degassing, dm; fuse=mdd.fuse)
        @. tracer += domain.step_tracer * mdd.volume_fraction[i]
    end
    # Cumulative fraction of tracer degassed
    cumsum!(fraction, tracer)
    fraction ./= last(fraction)
    return fraction
end
function degas_tracer!(sdd::SingleDomain{T}, dm::DiffusivityModel{T}) where {T<:AbstractFloat}
    fraction = fill!(sdd.model_fraction, zero(T))
    # Degas
    degas_tracer!(sdd.domain, one(T), sdd.tsteps_degassing, sdd.Tsteps_degassing, dm; fuse=sdd.fuse)
    # Cumulative fraction of tracer degassed
    cumsum!(fraction, sdd.domain.step_tracer)
    fraction ./= last(fraction)
    return fraction
end
function degas_tracer!(mineral::PlanarNobleGas{T}, initial_tracer, tsteps_degassing, Tsteps_degassing, dm; fuse::Bool=true)  where {T}
    # Erase previous diffusion profiles
    u = fill!(mineral.u, zero(T))
    u[2:end-1,1] .= initial_tracer
    u[1,1] = u[2,1] # Symmetric inner boundary condition
    u[end,1] = zero(T) # Zero concentration at outer boundary, since we degas in a vacuum
    # Degas and return total tracer
    diffusivityratio = tracerdiffusivityratio(mineral)
    crank_nicolson!(mineral, tsteps_degassing, Tsteps_degassing, dm; fuse, diffusivityratio)
    return diffusant_lost!(mineral.step_tracer, mineral; fuse)
end
function degas_tracer!(mineral::SphericalNobleGas{T}, initial_tracer, tsteps_degassing, Tsteps_degassing, dm; fuse::Bool=true)  where {T}
    # Erase previous diffusion profiles
    u = fill!(mineral.u, zero(T))
    u[2:end-1,1] .= initial_tracer
    u[2:end-1,1] .*= mineral.rsteps # U-transform for Crank-Nicholson
    u[1,1] = -u[2,1] # Symmetric inner boundary given U-transform
    u[end,1] = zero(T) # Zero concentration at outer boundary, since we degas in a vacuum
    # Degas and return total tracer
    diffusivityratio = tracerdiffusivityratio(mineral)
    crank_nicolson!(mineral, tsteps_degassing, Tsteps_degassing, dm; fuse, diffusivityratio)
    return diffusant_lost!(mineral.step_tracer, mineral; fuse)
end
tracerdiffusivityratio(x::ArgonSample{T}) where {T} = T((40/39)^0.3)
tracerdiffusivityratio(x::HeliumSample{T}) where {T} = T((4/3)^0.3)


## ---  Combined daughter+tracer degassing functions
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

function modelage(sdd::SingleDomain{T,<:HeliumSample}, Tsteps::AbstractVector, dm::DiffusivityModel{T}; redegastracer::Bool=false, partitiondaughter::Bool=false) where {T<:AbstractFloat}
    stepratio = fill!(sdd.model_age, zero(T))
    fraction = fill!(sdd.model_fraction, zero(T))
    # Degas
    age = modelage(sdd.domain, Tsteps, dm; partitiondaughter)
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
function modelage(sdd::SingleDomain{T,<:ArgonSample}, Tsteps::AbstractVector, dm::DiffusivityModel{T}; redegastracer::Bool=false, partitiondaughter::Bool=false) where {T<:AbstractFloat}
    stepage = fill!(sdd.model_age, zero(T))
    fraction = fill!(sdd.model_fraction, zero(T))
    # Degas
    age = modelage(sdd.domain, Tsteps, dm; partitiondaughter)
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
function modelage(mdd::MultipleDomain{T,<:ArgonSample}, Tsteps::AbstractVector, dm::MDDiffusivity{T}; redegastracer::Bool=false, partitiondaughter::Bool=false) where {T<:AbstractFloat}
    age = fill!(mdd.model_age, zero(T))
    tracer = fill!(mdd.model_tracer, zero(T))
    daughter = fill!(mdd.model_daughter, zero(T))
    fraction = fill!(mdd.model_fraction, zero(T))
    fuse = mdd.fuse::Bool
    # Degas
    for i in eachindex(mdd.domains, mdd.volume_fraction)
        domain = mdd.domains[i]
        modelage(domain, Tsteps, dm[i]; partitiondaughter)
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
        last_model_fractionᵢ = model_fractionᵢ
        last_fraction_experimentalᵢ = dd.fraction_experimental[i]
    end
    rescale && (ll /= sqrt(count(dd.fit)))
    return ll
end

## --- End of File
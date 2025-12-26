## --- Step heating diffusivity types

"""
```julia
SDDiffusivity(
    model::DiffusivityModel{T}  # Underlying (wrapped) diffusivity model
    scale::T                    # [unitless] relative domain size (default = 1.0)
    scale_logsigma::T           # [unitless] log uncertainty (default = 1.0 = a factor of ℯ, one-sigma)
)
```
One diffusivity, scaled to represent domain size as d/a^2
"""
Base.@kwdef struct SDDiffusivity{T,D<:DiffusivityModel{T}} <: DiffusivityModel{T}
    model::D                    # Underlying (wrapped) diffusivity model
    scale::T=1.0                # [unitless] relative domain size (default = 1.0)
    scale_logsigma::T=1.0       # [unitless] log uncertainty (default = 1.0 = a factor of ℯ, one-sigma)
end

"""
```julia
MDDiffusivity(
    D0::NTuple{N,T}             # [cm^2/sec] Maximum diffusivity
    D0_logsigma::NTuple{N,T}    # [unitless] log uncertainty
    Ea::T                       # [kJ/mol] Activation energy
    Ea_logsigma::T              # [unitless] log uncertainty
)
```
Multiple diffusivities for multiple domains
"""
Base.@kwdef struct MDDiffusivity{T<:AbstractFloat, N} <: DiffusivityModel{T}
    D0::NTuple{N,T}             # [cm^2/sec] Maximum diffusivity
    D0_logsigma::NTuple{N,T}    # [unitless] log uncertainty
    Ea::NTuple{N,T}             # [kJ/mol] Activation energy
    Ea_logsigma::NTuple{N,T}    # [unitless] log uncertainty 
end
Base.getindex(d::MDDiffusivity{T}, i::Int) where {T} = Diffusivity{T}(d.D0[i], d.D0_logsigma[i], d.Ea[i], d.Ea_logsigma[i])

# Implement eltype methods to deal with diffusivity models which are wrapper types
Base.eltype(x::DiffusivityModel) = typeof(x)
Base.eltype(x::SDDiffusivity{T,D}) where {T,D} = D

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

function degas!(mineral::HeliumSample, tsteps_degassing, Tsteps_degassing, dm; fuse::Bool=true, redegastracer::Bool=true)
    total_daughter = degas_daughter!(mineral, tsteps_degassing, Tsteps_degassing, dm; fuse)
    # Now diffuse parent isotope tracer (He-3), if neccesary
    if redegastracer || !(0 < sum(mineral.step_tracer))
        degas_tracer!(mineral, total_daughter, tsteps_degassing, Tsteps_degassing, dm; fuse)
    end
    return mineral.step_tracer, mineral.step_daughter
end
function degas!(mineral::ArgonSample, tsteps_degassing, Tsteps_degassing, dm; fuse::Bool=true, redegastracer::Bool=true)
    degas_daughter!(mineral, tsteps_degassing, Tsteps_degassing, dm; fuse)
    # Now diffuse parent isotope tracer (Ar-39), if neccesary
    if redegastracer || !(0 < sum(mineral.step_tracer))
        degas_tracer!(mineral, mineral.r40K, tsteps_degassing, Tsteps_degassing, dm; fuse)
    end
    return mineral.step_tracer, mineral.step_daughter
end


## --- Age and likelihood functions for step heating data

function modelage(sdd::SingleDomain{T,<:HeliumSample}, Tsteps::AbstractVector, dm::DiffusivityModel{T}; redegastracer::Bool=true, partitiondaughter::Bool=false) where {T<:AbstractFloat}
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
function modelage(sdd::SingleDomain{T,<:ArgonSample}, Tsteps::AbstractVector, dm::DiffusivityModel{T}; redegastracer::Bool=true, partitiondaughter::Bool=false) where {T<:AbstractFloat}
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
function modelage(mdd::MultipleDomain{T,<:ArgonSample}, Tsteps::AbstractVector, dm::MDDiffusivity{T}; redegastracer::Bool=true, partitiondaughter::Bool=false) where {T<:AbstractFloat}
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
function cumulative_fraction_uncertainty(sigma, i::Int)
    i₋ = firstindex(sigma)
    i₊ = lastindex(sigma)
    σ²₋ = sum(abs2, view(sigma, i₋:i))
    σ²₊ = sum(abs2, view(sigma,(i+1):i₊))
    σ = sqrt(1/(1/σ²₋ + 1/σ²₊))
end

## --- End of File
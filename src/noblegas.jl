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
# Apply to an ArgonSample
function newton_age(mineral::ArgonSample)
    # Daughter concentrations, in atoms/gram
    μAr = final_diffusant(mineral)
    # Parent concentrations in atoms/gram
    μ40K = meanparent(mineral)
    # Numerically solve for raw Ar age of the grain (i.e., as measured)
    return newton_ar_age(μAr, μ40K)
end
meanparent(mineral::PlanarAr) = nanmean(mineral.r40K)
meanparent(mineral::SphericalAr) = nanmean(mineral.r40K, mineral.relvolumes)

## --- Helium age functions

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
# Apply to a HeliumSample
function newton_age(mineral::HeliumSample)
    # Daughter concentrations, in atoms/gram
    μHe = final_diffusant(mineral)
    # Parent concentrations in atoms/gram
    μ238U, μ235U, μ232Th, μ147Sm = meanparent(mineral)
    # Numerically solve for raw Ar age of the grain (i.e., as measured)
    return newton_he_age(μHe, μ238U, μ235U, μ232Th, μ147Sm)
end
function meanparent(mineral::PlanarHe)
    return (nanmean(mineral.r238U), nanmean(mineral.r235U), nanmean(mineral.r232Th), nanmean(mineral.r147Sm))
end
function meanparent(mineral::Union{SphericalHe, ZirconHe, ApatiteHe})
    rv = mineral.relvolumes
    return (nanmean(mineral.r238U,rv), nanmean(mineral.r235U,rv), nanmean(mineral.r232Th,rv), nanmean(mineral.r147Sm,rv))
end

## --- Concrete types for damage and diffusivity models


"""
```julia
Diffusivity(
    D0::T = 59.98               # [cm^2/sec] Maximum diffusion coefficient
    D0_logsigma::T = log(2)/2   # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)
    Ea::T = 205.94              # [kJ/mol] Activation energy
    Ea_logsigma::T = log(2)/2   # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)
)
```
A generic diffusivity model, with user-specified D0 and Ea.
Default values are appropriate for argon in k-feldspar.
"""
Base.@kwdef struct Diffusivity{T<:AbstractFloat} <: DiffusivityModel{T}
    D0::T = 59.98               # [cm^2/sec] Maximum diffusion coefficient
    D0_logsigma::T = log(2)/2   # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)
    Ea::T = 205.94              # [kJ/mol] Activation energy
    Ea_logsigma::T = log(2)/2   # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)
end

"""
```julia
ZRDAAM(
    DzD0::T = 193188.0          # [cm^2/sec] Maximum diffusivity, crystalline endmember
    DzD0_logsigma::T=log(2)/2   # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)
    DzEa::T=165.0               # [kJ/mol] Activation energy, crystalline endmember
    DzEa_logsigma::T=log(2)/4   # [unitless] log uncertainty (default = log(2)/4 = a factor of sqrt(2) two-sigma)
    DN17D0::T = 6.367E-3        # [cm^2/sec] Maximum diffusivity, amorphous endmember
    DN17D0_logsigma::T=log(2)/2 # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)
    DN17Ea::T=71.0              # [kJ/mol] Activation energy, amorphous endmember
    DN17Ea_logsigma::T=log(2)/4 # [unitless] log uncertainty (default = log(2)/4 = a factor of sqrt(2) two-sigma)
    lint0::T=45920.0            # [nm]
    SV::T=1.669                 # [1/nm]
    Bα::T=5.48E-19              # Amorphous material produced per alpha decay [g/alpha]
    Phi::T=3.0                  # [unitless]
    beta::T=-0.05721            # Zircon anealing parameter
    C0::T=6.24534               # Zircon anealing parameter
    C1::T=-0.11977              # Zircon anealing parameter
    C2::T=-314.937 - LOG_SEC_MYR # Zircon anealing parameter. Includes conversion factor from seconds to Myr for dt (for performance), in addition to traditional C2 value
    C3::T=-14.2868              # Zircon anealing parameter
    rmin::T=0.2                 # Damage conversion parameter
    rmin_sigma::T=0.15          # Damage conversion parameter uncertainty
)
```
Zircon Radiation Damage Accumulation and Annealing Model (ZRDAAM) of
Guenthner et al. 2013 (doi: 10.2475/03.2013.01)
"""
Base.@kwdef struct ZRDAAM{T<:AbstractFloat} <: ZirconHeliumModel{T} 
    DzD0::T = 193188.0          # [cm^2/sec] Maximum diffusivity, crystalline endmember
    DzD0_logsigma::T=log(2)/2   # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)
    DzEa::T=165.0               # [kJ/mol] Activation energy, crystalline endmember
    DzEa_logsigma::T=log(2)/4   # [unitless] log uncertainty (default = log(2)/4 = a factor of sqrt(2) two-sigma)
    DN17D0::T = 6.367E-3        # [cm^2/sec] Maximum diffusivity, amorphous endmember
    DN17D0_logsigma::T=log(2)/2 # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)
    DN17Ea::T=71.0              # [kJ/mol] Activation energy, amorphous endmember
    DN17Ea_logsigma::T=log(2)/4 # [unitless] log uncertainty (default = log(2)/4 = a factor of sqrt(2) two-sigma)
    lint0::T=45920.0            # [nm]
    SV::T=1.669                 # [1/nm]
    Bα::T=5.48E-19              # Amorphous material produced per alpha decay [g/alpha]
    Phi::T=3.0                  # [unitless]
    beta::T=-0.05721            # Zircon anealing parameter
    C0::T=6.24534               # Zircon anealing parameter
    C1::T=-0.11977              # Zircon anealing parameter
    C2::T=-314.937 - LOG_SEC_MYR # Zircon anealing parameter. Includes conversion factor from seconds to Myr for dt (for performance), in addition to traditional C2 value
    C3::T=-14.2868              # Zircon anealing parameter
    rmin::T=0.2                 # Damage conversion parameter
    rmin_sigma::T=0.15          # Damage conversion parameter uncertainty
end

"""
```julia
RDAAM(
    D0L::T=0.6071               # [cm^2/sec] Maximum diffusivity
    D0L_logsigma::T=log(2)/2    # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)
    EaL::T=122.3                # [kJ/mol] Activation energy
    EaL_logsigma::T=log(2)/4    # [unitless] log uncertainty (default = log(2)/4 = a factor of sqrt(2) two-sigma)
    EaTrap::T=34.0              # [kJ/mol] Activation energy
    EaTrap_logsigma::T=log(2)/4 # [unitless] log uncertainty (default = log(2)/4 = a factor of sqrt(2) two-sigma)
    psi::T=1e-13                # empirical polynomial coefficient
    omega::T=1e-22              # empirical polynomial coefficient
    etaq::T=0.91                # Durango ηq
    rhoap::T=3.19               # [g/cm3] Density of apatite 
    L::T=0.000815               # [cm] Etchable fission track half-length
    lambdaf::T=8.46e-17         # [1/a] U-238 spontaneous fission decay constant
    lambdaD::T=1.55125e-10      # [1/a] U-238 total decay constant
    beta::T=0.04672             # Apatite annealing parameter. Also caled alpha, but equivalent to beta in ZRDAAM
    C0::T=0.39528               # Apatite annealing parameter
    C1::T=0.01073               # Apatite annealing parameter
    C2::T=-65.12969 - LOG_SEC_MYR # Apatite annealing parameter. Includes conversion factor from seconds to Myr for dt, in addition to traditional C2 value
    C3::T=-7.91715              # Apatite annealing parameter
    rmr0::T=0.83                # Damage conversion parameter
    rmr0_sigma::T=0.15          # Damage conversion parameter uncertainty
    kappa::T=1.04-0.83          # Damage conversion parameter
    kappa_rmr0::T=1.04          # Damage conversion parameter (the sum of kappa and rmr0)
)
```
Apatite Radiation Damage Accumulation and Annealing Model (RDAAM) of
Flowers et al. 2009 (doi: 10.1016/j.gca.2009.01.015)
"""
Base.@kwdef struct RDAAM{T<:AbstractFloat} <: ApatiteHeliumModel{T} 
    D0L::T=0.6071               # [cm^2/sec] Maximum diffusivity
    D0L_logsigma::T=log(2)/2    # [unitless] log uncertainty (default = log(2)/2 = a factor of 2 two-sigma)
    EaL::T=122.3                # [kJ/mol] Activation energy
    EaL_logsigma::T=log(2)/4    # [unitless] log uncertainty (default = log(2)/4 = a factor of sqrt(2) two-sigma)
    EaTrap::T=34.0              # [kJ/mol] Activation energy
    EaTrap_logsigma::T=log(2)/4 # [unitless] log uncertainty (default = log(2)/4 = a factor of sqrt(2) two-sigma)    psi::T=1e-13                # empirical polynomial coefficient
    psi::T=1e-13                # empirical polynomial coefficient
    omega::T=1e-22              # empirical polynomial coefficient
    etaq::T=0.91                # Durango ηq
    rhoap::T=3.19               # [g/cm3] Density of apatite 
    L::T=0.000815               # [cm] Etchable fission track half-length
    lambdaf::T=8.46e-17         # [1/a] U-238 spontaneous fission decay constant
    lambdaD::T=1.55125e-10      # [1/a] U-238 total decay constant
    beta::T=0.04672             # Apatite annealing parameter. Also caled alpha, but equivalent to beta in ZRDAAM
    C0::T=0.39528               # Apatite annealing parameter
    C1::T=0.01073               # Apatite annealing parameter
    C2::T=-65.12969 - LOG_SEC_MYR # Apatite annealing parameter. Includes conversion factor from seconds to Myr for dt (for performance), in addition to traditional C2 value
    C3::T=-7.91715              # Apatite annealing parameter
    rmr0::T=0.83                # Damage conversion parameter
    rmr0_sigma::T=0.15          # Damage conversion parameter uncertainty
    kappa::T=1.04-0.83          # Damage conversion parameter
    kappa_rmr0::T=1.04          # Damage conversion parameter (the sum of kappa and rmr0)
end

# Functions for querying helium models from kinetic results
function ZirconHeliumModel(kr::KineticResult)
    any(x->isa(x, ZirconHeliumModel), kr) || return nothing
    i = findfirst(x->isa(x, ZirconHeliumModel), kr[:,1])
    return vec(kr[i,:])
end
function ApatiteHeliumModel(kr::KineticResult)
    any(x->isa(x, ApatiteHeliumModel), kr) || return nothing
    i = findfirst(x->isa(x, ApatiteHeliumModel), kr[:,1])
    return vec(kr[i,:])
end

## --- Damage and annealing functions

"""
```julia
ρᵣ = anneal(dt::Number, tsteps::Vector, Tsteps::Matrix, [model::DiffusivityModel=ZRDAAM()])
```
ZirconHe damage annealing model as in Guenthner et al. 2013 (AJS)
"""
function anneal(tsteps::AbstractVector, Tsteps::AbstractVector, dm::DiffusivityModel=ZRDAAM())
    # Allocate matrix to hold reduced track lengths for all previous timesteps
    ntsteps = length(tsteps)
    ρᵣ = zeros(ntsteps,ntsteps)
    teq = zeros(ntsteps)
    # In=-place version
    anneal!(ρᵣ, teq, tsteps, Tsteps, dm)
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
function anneal!(data::Vector{<:Chronometer{T}}, ::Type{C}, tsteps::AbstractVector{T}, Tsteps::AbstractVector{T}, dm::DiffusivityModel{T}) where {T<:AbstractFloat, C<:HeliumSample}
    @assert eachindex(tsteps) == eachindex(Tsteps)
    if any(x->isa(x, C), data)
        im = argmax(i->(eltype(data[i]) <: C) ? length(timediscretization(data[i])) : 0, eachindex(data))
        first_index = firstindex(Tsteps) + length(tsteps) - length(timediscretization(data[im]))
        dₘ = if first_index > 1
            anneal!(data[im], @views(Tsteps[first_index:end]), dm)::C
        else
            anneal!(data[im], Tsteps, dm)::C
        end
        pr = dₘ.pr
        @assert length(timediscretization(dₘ)) == length(axes(pr, 1)) == length(axes(pr, 2))
        for i in eachindex(data)
            if i!=im && eltype(data[i]) <: C
                anneal!(data[i], Tsteps, pr)
            end
        end
    end
    return data
end
function anneal!(mineral::HeliumSample, Tsteps::AbstractVector, pr::AbstractMatrix)
    ntsteps = length(axes(pr, 1))
    first_index = firstindex(pr) + ntsteps - length(timediscretization(mineral))
    if first_index > 1
        mul!(mineral.annealeddamage, @views(pr[first_index:end, first_index:end]), mineral.alphadamage)
    else
        mul!(mineral.annealeddamage, pr, mineral.alphadamage)
    end
    return mineral
end
anneal!(mineral::SingleDomain, Tsteps::AbstractVector, prdm) = anneal!(mineral.domain, Tsteps, prdm)
function anneal!(mineral::MultipleDomain, Tsteps::AbstractVector, prdm)
    c₀ = first(mineral.domains)
    anneal!(c₀, Tsteps, prdm)
    pr = c₀.pr
    for c in Iterators.drop(mineral.domains,1)
        mul!(c.annealeddamage, pr, c.alphadamage)
    end
    return c₀
end
function anneal!(mineral::ZirconHe, Tsteps::AbstractVector, dm::ZRDAAM)
    anneal!(mineral.pr, view(mineral.annealeddamage,:,1), mineral.tsteps, Tsteps, dm)
    mul!(mineral.annealeddamage, mineral.pr, mineral.alphadamage)
    return mineral
end
function anneal!(mineral::ApatiteHe, Tsteps::AbstractVector, dm::RDAAM)
    anneal!(mineral.pr, view(mineral.annealeddamage,:,1), mineral.tsteps, Tsteps, dm)
    mul!(mineral.annealeddamage, mineral.pr, mineral.alphadamage)
    return mineral
end
function anneal!(ρᵣ::AbstractMatrix{T}, teq::AbstractVector{T}, tsteps::AbstractVector, Tsteps::AbstractVector, dm::ZRDAAM{T}) where T <: AbstractFloat
    @assert eachindex(tsteps) == eachindex(Tsteps) == eachindex(teq) == axes(ρᵣ, 1) == axes(ρᵣ,2) == Base.OneTo(length(tsteps))
    ntsteps = length(tsteps)
    fill!(teq, zero(T))

    # First timestep
    dt = step_at(tsteps, 1)
    ρᵣ[1,1] = 1 / ((dm.C0 + dm.C1*(log(dt)-dm.C2)/(log(1 / (Tsteps[1]+273.15))-dm.C3))^(1/dm.beta)+1)

    # All subsequent timesteps
    @inbounds for i=2:ntsteps
        dt = step_at(tsteps, i)
        lᵢ = log(1 / (Tsteps[i]+273.15)) - dm.C3

        # Convert any existing track length reduction for damage from
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
function anneal!(ρᵣ::AbstractMatrix{T}, teq::AbstractVector{T}, tsteps::AbstractVector, Tsteps::AbstractVector, dm::RDAAM{T}) where T <: AbstractFloat

    ntsteps = length(tsteps)
    @assert size(ρᵣ) === (ntsteps, ntsteps)
    @assert size(teq) === (ntsteps,)
    fill!(teq, zero(T))

    # First timestep
    dt = step_at(tsteps, 1)
    ρᵣ[1,1] = 1 / ((dm.C0 + dm.C1*(log(dt)-dm.C2)/(log(1 / (Tsteps[1]+273.15))-dm.C3))^(1/dm.beta)+1)

    # All subsequent timesteps
    @inbounds for i=2:ntsteps
        dt = step_at(tsteps, i)
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

## --- Age and likelihood of a noble gas chronometer given a t-T path

"""
```julia
modelage(mineral::ZirconHe, Tsteps, [ρᵣ], dm::ZRDAAM)
modelage(mineral::ApatiteHe, Tsteps, [ρᵣ], dm::RDAAM)
modelage(mineral::SphericalHe, Tsteps, dm::Diffusivity)
modelage(mineral::PlanarHe, Tsteps, dm::Diffusivity)
modelage(mineral::SphericalAr, Tsteps, dm::Diffusivity)
modelage(mineral::PlanarAr, Tsteps, dm::Diffusivity)
```
Calculate the predicted bulk age of a noble gas chronometer that has experienced a given 
t-T path (specified by `mineral.tsteps` for time and `Tsteps` for temperature), 
at a time resolution determined by `mineral.tsteps` using a Crank-Nicolson diffusion 
solution for a spherical (or planar slab) grain of radius (or halfwidth ) `mineral.r` 
at spatial resolution `mineral.dr`.

Spherical implementation based on the the Crank-Nicolson solution for diffusion out of a
spherical mineral crystal in Ketcham, 2005 (doi: 10.2138/rmg.2005.58.11).
"""
function modelage(mineral::NobleGasSample{T}, Tsteps::AbstractVector, dm::DiffusivityModel{T}; partitiondaughter::Bool=false) where {T <: AbstractFloat}
    # Erase any previous runs; start at zero initial daughter
    fill!(mineral.u, zero(T))

    # Run Crank-Nicolson solver
    crank_nicolson!(mineral, mineral.tsteps, Tsteps, dm; partitiondaughter, fuse=false, setting=:geological) 

    # Numerically solve for resulting observed age of the grain (i.e, as measured, "raw" in AHe/ZHe parlanc)
    return newton_age(mineral)::T
end
function modelage(mineral::Union{ZirconHe,ApatiteHe}, Tsteps::AbstractVector, ρᵣ::AbstractMatrix, dm::Union{RDAAM,ZRDAAM})
    mul!(mineral.annealeddamage, ρᵣ, mineral.alphadamage)
    modelage(mineral, Tsteps, dm)
end

# Log likelihood for model ages
function model_ll(mineral::NobleGasSample, Tsteps, dm)
    age = modelage(mineral, Tsteps, dm)
    δ = age - mineral.age
    σ² = mineral.age_sigma^2
    -0.5*(log(2*pi*σ²) + δ^2/σ²)
end
function model_ll(mineral::Union{ZirconHe,ApatiteHe}, Tsteps, dm::DiffusivityModel)
    anneal!(mineral, Tsteps, dm)
    age = modelage(mineral, Tsteps, dm)
    δ = age - mineral.age
    σ² = mineral.age_sigma^2
    -0.5*(log(2*pi*σ²) + δ^2/σ²)
end

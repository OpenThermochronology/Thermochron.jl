## --- Define AnnealingModel types

abstract type AnnealingModel{T} end

# Implement methods to allow broadcasting
Base.length(x::AnnealingModel) = 1
Base.iterate(x::AnnealingModel) = (x, nothing)
Base.iterate(x::AnnealingModel, state) = nothing

struct FanningCurvilinear{T<:AbstractFloat} <: AnnealingModel{T} 
    C0::T
    C1::T
    C2::T
    C3::T
    alpha::T
    beta::T
    l0::T
    l0_sigma::T
end

struct SimplifiedCurvilinear{T<:AbstractFloat} <: AnnealingModel{T} 
    C0::T
    C1::T
    C2::T
    C3::T
    alpha::T
    l0::T
    l0_sigma::T
end

## --- Define DiffusivityModel types

abstract type DiffusivityModel{T} end
abstract type ZirconHeliumModel{T} <: DiffusivityModel{T} end
abstract type ApatiteHeliumModel{T} <: DiffusivityModel{T} end

# Implement methods to allow broadcasting
Base.length(x::DiffusivityModel) = 1
Base.iterate(x::DiffusivityModel) = (x, nothing)
Base.iterate(x::DiffusivityModel, state) = nothing

Base.@kwdef struct ZRDAAM{T<:AbstractFloat} <: ZirconHeliumModel{T} 
    DzD0::T = 193188.0          # Diffusivity [cm^2/sec], crystalline endmember
    DzD0_logsigma::T=log(2)     # log units (default = log(2) = a factor of 2)
    DzEa::T=165.0               # Activation energy [kJ/mol], crystalline endmember
    DzEa_logsigma::T=log(2)     # log units (default = log(2) = a factor of 2)
    DN17D0::T = 6.367E-3        # Diffusivity [cm^2/sec], amorphous endmember
    DN17D0_logsigma::T=log(2)   # log units (default = log(2) = a factor of 2)
    DN17Ea::T=71.0              # Activation energy [kJ/mol], amorphous endmember
    DN17Ea_logsigma::T=log(2)   # log units (default = log(2) = a factor of 2)
    lint0::T=45920.0            # [nm]
    SV::T=1.669                 # [1/nm]
    Bα::T=5.48E-19              # Amorphous material produced per alpha decay [g/alpha]
    Phi::T=3.0                  # unitless
    beta::T=-0.05721            # Zircon anealing parameter
    C0::T=6.24534               # Zircon anealing parameter
    C1::T=-0.11977              # Zircon anealing parameter
    C2::T=-314.937 - LOG_SEC_MYR # Zircon anealing parameter. Includes conversion factor from seconds to Myr for dt (for performance), in addition to traditional C2 value
    C3::T=-14.2868              # Zircon anealing parameter
    rmr0::T=0.2                 # Damage conversion parameter (not normally called this in ZRDAAM context, but effectively the same thing)
    rmr0_sigma::T=0.1           # Damage conversion parameter uncertainty
end
export ZRDAAM

Base.@kwdef struct RDAAM{T<:AbstractFloat} <: ApatiteHeliumModel{T} 
    D0L::T=0.6071               # Diffusivity [cm^2/s]
    D0L_logsigma::T=log(2)      # log units (default = log(2) = a factor of 2)
    EaL::T=122.3                # Activation energy [kJ/mol]
    EaL_logsigma::T=log(2)      # log units (default = log(2) = a factor of 2)
    EaTrap::T=34.0              # Activation energy [kJ/mol]
    EaTrap_logsigma::T=log(2)   # log units (default = log(2) = a factor of 2)
    psi::T=1e-13                # empirical polynomial coefficient
    omega::T=1e-22              # empirical polynomial coefficient
    etaq::T=0.91                # Durango ηq
    rhoap::T=3.19               # Density of apatite [g/cm3]
    L::T=0.000815               # Etchable fission track half-length, cm
    lambdaf::T=8.46e-17         # 
    lambdaD::T=1.55125e-10      # 
    beta::T=0.04672             # Apatite annealing parameter. Also caled alpha, but equivalent to beta in ZRDAAM
    C0::T=0.39528               # Apatite annealing parameter
    C1::T=0.01073               # Apatite annealing parameter
    C2::T=-65.12969 - LOG_SEC_MYR # Apatite annealing parameter. Includes conversion factor from seconds to Myr for dt (for performance), in addition to traditional C2 value
    C3::T=-7.91715              # Apatite annealing parameter
    rmr0::T=0.83                # Damage conversion parameter
    rmr0_sigma::T=0.1           # Damage conversion parameter uncertainty
    kappa::T=1.04-0.83          # Damage conversion parameter
    kappa_rmr0::T=1.04          # Damage conversion parameter (the sum of kappa and rmr0)
end
export RDAAM

## --- Define Boundary type to specify the working area

struct Boundary{T<:AbstractFloat}
    agepoints::NTuple{2,T}  # Ma
    Tpoints::Vector{T}      # Degrees C
    Tpointsₚ::Vector{T}     # Degrees C
    T₀::NTuple{2,T}         # Degrees C
    ΔT::NTuple{2,T}         # Degrees C
    tboundary::Symbol
    Tboundary::Symbol
    npoints::Int
end
function Boundary(T::Type=Float64; agepoints, T₀, ΔT, Tpoints=collect(T₀), tboundary=:reflecting, Tboundary=:reflecting)
    @assert length(agepoints) == length(T₀) == length(ΔT) == length(Tpoints) == 2
    Boundary{T}(
        Tuple(T.(agepoints)), 
        collect(T, Tpoints), 
        collect(T, Tpoints), 
        Tuple(T.(T₀)), 
        Tuple(T.(ΔT)), 
        Symbol(tboundary),
        Symbol(Tboundary),
        2,
    )
end
export Boundary

## -- Define Constraint type to specify aribitrary t-T constraint boxes
struct Constraint{T<:AbstractFloat}
    agepoints::Vector{T} # Ma
    Tpoints::Vector{T}   # Degrees C
    agepointsₚ::Vector{T} # Ma
    Tpointsₚ::Vector{T}   # Degrees C
    agedist::Vector{<:Union{Uniform{T}, Normal{T}, LogNormal{T}}}
    Tdist::Vector{<:Union{Uniform{T}, Normal{T}, LogNormal{T}}}
    npoints::Int
end
function Constraint(T::Type=Float64; agedist=Uniform{T}[], Tdist=Uniform{T}[], agepoints=mean.(agedist), Tpoints=mean.(Tdist))
    @assert length(agepoints) == length(Tpoints) == length(agedist) == length(Tdist)
    isabstracttype(eltype(agedist)) && (agedist = unionize(agedist))
    isabstracttype(eltype(Tdist)) && (Tdist = unionize(Tdist))
    Constraint{T}(T.(agepoints),
        T.(Tpoints),
        T.(agepoints), 
        T.(Tpoints), 
        agedist,
        Tdist,
        length(agepoints),
    )
end
export Constraint

# For backwards compatibility with old scripts
const Unconformity = Constraint
export Unconformity

## --- Define DetailInterval type to specify a minumum number of t-T path nodes within a given time interval
struct DetailInterval{T<:AbstractFloat}
    agemin::T
    agemax::T
    minpoints::Int
end
function DetailInterval(T::Type=Float64; agemin=0, agemax=0, minpoints=0)
    DetailInterval{T}(agemin, agemax, Int(minpoints))
end
export DetailInterval


## --- Model result types

abstract type AbstractTTResult end

struct TTResult{T<:AbstractFloat} <: AbstractTTResult
    tpointdist::Matrix{T}
    Tpointdist::Matrix{T}
    ndist::Vector{Int}
    HeAgedist::Matrix{T}
    jtdist::Vector{T}
    jTdist::Vector{T}
    lldist::Vector{T}
    acceptancedist::BitVector
end

abstract type AbstractKineticResult end

struct KineticResult{T<:AbstractFloat} <: AbstractKineticResult
    admdist::Vector{<:ApatiteHeliumModel{T}}
    zdmdist::Vector{<:ZirconHeliumModel{T}}
end
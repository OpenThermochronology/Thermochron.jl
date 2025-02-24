## --- Define AnnealingModel type hierarchy

abstract type Model{T} end
abstract type AnnealingModel{T} <: Model{T} end
abstract type ZirconAnnealingModel{T} <: AnnealingModel{T} end
abstract type MonaziteAnnealingModel{T} <: AnnealingModel{T} end
abstract type ApatiteAnnealingModel{T} <: AnnealingModel{T} end
abstract type FanningCurvilinearZircon{T} <: ZirconAnnealingModel{T} end
abstract type FanningCurvilinearApatite{T} <: ApatiteAnnealingModel{T} end
const FanningCurvilinear{T} = Union{FanningCurvilinearZircon{T}, FanningCurvilinearApatite{T}}

# Implement methods to allow broadcasting and comparison
Base.length(x::Model) = 1
Base.iterate(x::Model) = (x, nothing)
Base.iterate(x::Model, state) = nothing
Base.:(==)(x::Model, y::Model) = false
function Base.:(==)(x::T, y::T) where {T<:Model}
    for n in fieldnames(T)
        isequal(getfield(x, n), getfield(y, n)) || return false
    end
    return true
end

## --- Define DiffusivityModel type hierarchy

abstract type DiffusivityModel{T} <: Model{T} end
abstract type ZirconHeliumModel{T} <: DiffusivityModel{T} end
abstract type ApatiteHeliumModel{T} <: DiffusivityModel{T} end

# Implement methods to allow broadcasting
Base.length(x::DiffusivityModel) = 1
Base.iterate(x::DiffusivityModel) = (x, nothing)
Base.iterate(x::DiffusivityModel, state) = nothing


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

# For backwards compatibility with old scripts
const Unconformity = Constraint

# Define DetailInterval type to specify a minumum number of t-T path nodes within a given time interval
struct DetailInterval{T<:AbstractFloat}
    agemin::T
    agemax::T
    minpoints::Int
end
function DetailInterval(T::Type=Float64; agemin=0, agemax=0, minpoints=0)
    DetailInterval{T}(agemin, agemax, Int(minpoints))
end

# Define overall TtPath type to contain all parameters needed to construct a t-T path proposal
struct TtPath{T<:AbstractFloat}
    agesteps::FloatRange
    Tsteps::Vector{T}
    agepoints::Vector{T}
    Tpoints::Vector{T}
    agepointsₚ::Vector{T}
    Tpointsₚ::Vector{T}
    σⱼt::Vector{T}
    σⱼT::Vector{T}
    σⱼtₚ::Vector{T}
    σⱼTₚ::Vector{T}
    agepointbuffer::Vector{T}
    Tpointbuffer::Vector{T}
    knot_index::Vector{Int}
    constraint::Constraint{T}
    boundary::Boundary{T}
    detail::DetailInterval{T}
end
function TtPath(agesteps::AbstractArray, constraint::Constraint{T}, boundary::Boundary{T}, detail::DetailInterval{T}, maxpoints::Int) where {T}
    # Discretized temperature
    agesteps = floatrange(agesteps)
    Tsteps = zeros(T, length(agesteps))
    knot_index = zeros(Int, length(agesteps))

    # Arrays to hold all t and T points (up to npoints=maxpoints)
    agepoints = zeros(T, maxpoints) 
    Tpoints = zeros(T, maxpoints)
    agepointsₚ = zeros(T, maxpoints) 
    Tpointsₚ = zeros(T, maxpoints)

    # Arrays to hold standard deviations of Gaussian proposal ("jumping") distributions for t and T
    σⱼt = fill(nanrange(textrema(boundary))/60, maxpoints)
    σⱼT = fill(nanrange(Textrema(boundary))/60, maxpoints)
    σⱼtₚ = copy(σⱼt)
    σⱼTₚ = copy(σⱼT)

    # Calculate number of boundary and unconformity points and allocate pointbuffer for interpolating
    totalpoints = maxpoints + boundary.npoints + constraint.npoints
    agepointbuffer = similar(agepoints, totalpoints)
    Tpointbuffer = similar(agepoints, totalpoints)
    TtPath{T}(
        agesteps,
        Tsteps,
        agepoints,
        Tpoints,
        agepointsₚ,
        Tpointsₚ,
        σⱼt,
        σⱼT,
        σⱼtₚ,
        σⱼTₚ,
        agepointbuffer,
        Tpointbuffer,
        knot_index,
        constraint,
        boundary,
        detail,
    )
end

## --- Model result types

abstract type AbstractTTResult end

struct TTResult{T<:AbstractFloat} <: AbstractTTResult
    tpointdist::Matrix{T}
    Tpointdist::Matrix{T}
    ndist::Vector{Int}
    resultdist::Matrix{T}
    jtdist::Vector{T}
    jTdist::Vector{T}
    lldist::Vector{T}
    acceptancedist::BitVector
end

abstract type AbstractKineticResult end

struct KineticResult{T<:AbstractFloat} <: AbstractKineticResult
    dmdist::Matrix{<:Model{T}}
end

function zdmdist(kr::KineticResult)
    any(x->isa(x, ZirconHeliumModel), kr) || return nothing
    i = findfirst(x->isa(x, ZirconHeliumModel), kr[:,1])
    return vec(kr[i,:])
end
function admdist(kr::KineticResult)
    any(x->isa(x, ApatiteHeliumModel), kr) || return nothing
    i = findfirst(x->isa(x, ApatiteHeliumModel), kr[:,1])
    return vec(kr[i,:])
end
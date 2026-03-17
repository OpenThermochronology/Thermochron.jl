## --- Chronometer type hierarchy

# Abstract type to include any number of mineral chronometers (zircon, apatite, etc.)
abstract type Chronometer{T} end

# Implement methods to allow broadcasting
Base.length(x::Chronometer) = 1
Base.iterate(x::Chronometer) = (x, nothing)
Base.iterate(x::Chronometer, state) = nothing

# Implement methods to allow copying and comparison
Base.copy(x::Chronometer) = deepcopy(x)
Base.:(==)(x::Chronometer, y::Chronometer) = false
function Base.:(==)(x::T, y::T) where {T<:Chronometer}
    for n in fieldnames(T)
        isequal(getfield(x, n), getfield(y, n)) || return false
    end
    return true
end

# Abstract subtypes for Chronometers that include an absolute age and age uncertainty
abstract type AbsoluteChronometer{T} <:Chronometer{T} end  
abstract type FissionTrackSample{T} <: AbsoluteChronometer{T} end   # Any absolute chronometer based on the annealing of fission tracks
abstract type StepHeatingSample{T,C} <: AbsoluteChronometer{T} end  # A chronometer recording laboratory step heating. Wraps one or more diffusive (typically noble gas) chronometers
abstract type NobleGasSample{T} <: AbsoluteChronometer{T} end       # Any chronometer based on diffusion of noble gasses
abstract type HeliumSample{T} <: NobleGasSample{T} end              # Any chronometer based on the diffusion of radiogenic He
abstract type ArgonSample{T} <: NobleGasSample{T} end               # Any chronometer based on the diffusion of radiogenic Ar

# Other abstract Chronometer types
abstract type FissionTrackLength{T} <: Chronometer{T} end           # Any chronometer based on the lengths of partially annealed fission tracks

## --- Model type hierarchy (annealing and diffusion models)
abstract type Model{T} end

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
Base.isapprox(x::Model, y::Model) = false
function Base.isapprox(x::T, y::T; kwargs...) where {T<:Model}
    for n in fieldnames(T)
        isapprox(getfield(x, n), getfield(y, n), kwargs...) || return false
    end
    return true
end

# Methods to allow averaging of models
Base.zero(::M) where {M<:Model} = zero(M)
@generated function Base.zero(::Type{M}) where {T, M<:Model{T}}
    result = :($M())
    for e in fieldnames(M)
        push!(result.args, :(zero(T)))
    end
    return result
end
@generated function Base.:+(x::M, y::M) where {M<:Model}
    result = :($M())
    for e in fieldnames(M)
        push!(result.args, :(x.$e + y.$e))
    end
    return result
end
@generated function Base.:*(x::M, n::Number) where {M<:Model}
    result = :($M())
    for e in fieldnames(M)
        push!(result.args, :(x.$e * n))
    end
    return result
end
Base.:*(n::Number, x::Model) = x * n          # Scalar multiplication is commutative
Base.:/(x::Model, n::Number) = x * inv(n)     # Division by a scalar is multiplciation by multiplicative inverse


# AnnealingModel types
abstract type AnnealingModel{T} <: Model{T} end
abstract type ZirconAnnealingModel{T} <: AnnealingModel{T} end
abstract type MonaziteAnnealingModel{T} <: AnnealingModel{T} end
abstract type ApatiteAnnealingModel{T} <: AnnealingModel{T} end
abstract type FanningCurvilinearZircon{T} <: ZirconAnnealingModel{T} end
abstract type FanningCurvilinearApatite{T} <: ApatiteAnnealingModel{T} end
const FanningCurvilinear{T} = Union{FanningCurvilinearZircon{T}, FanningCurvilinearApatite{T}}

# DiffusivityModel types
abstract type DiffusivityModel{T} <: Model{T} end
abstract type MultipleDiffusivity{T} <: DiffusivityModel{T} end
abstract type ZirconHeliumModel{T} <: DiffusivityModel{T} end
abstract type ApatiteHeliumModel{T} <: DiffusivityModel{T} end

# Other model types
Base.@kwdef struct RegionalParameters{T} <: Model{T}
    geotherm::T = 30.0                  # [C/km] Geothermal gradient
    geotherm_logsigma::T = log(2)       # [unitless] Log uncertainty, factor of 2 (one-sigma)
    K0_itm::T = 0.00010632              # [unitless] Preexponential constant for mineral-ITM/IGB noble gas partitioning
    K0_itm_logsigma::T = log(2)         # [unitless] log uncertainty, factor of 2 (one-sigma)
    Ea_itm::T = 26.815                  # [kJ/mol] Activation energy for for mineral-ITM/IGB noble gas partitioning
    Ea_itm_logsigma::T = log(2)/2       # [unitless] log uncertainty, factor of sqrt(2) (one-sigma)
    Ea_lambda::T = 30.0                 # [kJ/mol] Effective activation energy associated with bulk loss rate λ [1/Ma]
    Ea_lambda_logsigma::T = log(2)/2    # [unitless] log uncertainty, factor of sqrt(2) (one-sigma)
end

## --- Types used in t-T path generation

# Define Boundary types to specify the working area
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

# Define Constraint type to specify aribitrary t-T constraint boxes
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
function Unconformity(args...)
    @warn "The `Unconformity` type has been deprecated, use `Constraint` instead"
    Constraint(args...)
end
export Unconformity

# Define DetailInterval type to specify a minumum number of t-T path nodes within a given time interval
struct DetailInterval{T<:AbstractFloat}
    agemin::T
    agemax::T
    minpoints::Int
end
function DetailInterval(T::Type=Float64; agemin=0, agemax=0, minpoints=0)
    DetailInterval{T}(agemin, agemax, Int(minpoints))
end

# Define overall TTPath type to contain all parameters needed to construct a t-T path proposal
struct TTPath{T<:AbstractFloat, V<:AbstractVector{T}}
    agesteps::V
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
function TTPath(agesteps::AbstractVector{T}, constraint::Constraint{T}, boundary::Boundary{T}, detail::DetailInterval{T}, maxpoints::Int) where {T}
    # Discretized temperature
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
    TTPath(
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

struct TTResult{T<:AbstractFloat, V<:AbstractVector{T}} <: AbstractVector{Vector{T}}
    agesteps::V
    tpointdist::Matrix{T}
    Tpointdist::Matrix{T}
    ndist::Vector{Int}
    resultdist::Matrix{T}
    jtdist::Vector{T}
    jTdist::Vector{T}
    lldist::Vector{T}
    acceptancedist::BitVector
end
Base.size(tT::TTResult) = (size(tT.tpointdist, 2),)
function Base.getindex(tT::TTResult{T}, i::Int) where {T}
    ages = view(tT.tpointdist, 1:(tT.ndist[i]+2))
    temperatures = view(tT.Tpointdist, 1:(tT.ndist[i]+2))
    return linterp1s(ages, temperatures, tT.agesteps) # Tsteps
end

struct KineticResult{T<:AbstractFloat, M<:Model{T}} <: AbstractVector{SubArray{M, 1, Matrix{M}, Tuple{Base.Slice{Base.OneTo{Int}}, Int}, true}}
    dmdist::Matrix{M}
end
Base.size(k::KineticResult) = (size(k.dmdist, 2),)
Base.getindex(k::KineticResult, i::Int) = view(k.dmdist, :, i)

## --- End of File
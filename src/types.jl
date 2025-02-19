## --- Define AnnealingModel types

abstract type AnnealingModel{T} end
abstract type ZirconAnnealingModel{T} <: AnnealingModel{T} end
abstract type MonaziteAnnealingModel{T} <: AnnealingModel{T} end
abstract type ApatiteAnnealingModel{T} <: AnnealingModel{T} end
abstract type FanningCurvilinearZircon{T} <: ZirconAnnealingModel{T} end
abstract type FanningCurvilinearApatite{T} <: ApatiteAnnealingModel{T} end
const FanningCurvilinear{T} = Union{FanningCurvilinearZircon{T}, FanningCurvilinearApatite{T}}

# Implement methods to allow broadcasting
Base.length(x::AnnealingModel) = 1
Base.iterate(x::AnnealingModel) = (x, nothing)
Base.iterate(x::AnnealingModel, state) = nothing

"""
```julia
Yamada2007PC(
    c0p = -63.37     # Yamada et al. 2007 zircon
    c1p = 0.212      # Yamada et al. 2007 zircon
    bp = 43.00       # Yamada et al. 2007 zircon
    l0 = 11.17       # [um] effective initial track length (μmax)
    l0_sigma = 0.051 # [um] effective initial track length uncertainty (σ)
)
```
Parallel Curvilinear zircon annealing model of Yamada, 2007 
(doi: 10.1016/j.chemgeo.2006.09.002)
"""
Base.@kwdef struct Yamada2007PC{T<:AbstractFloat} <: ZirconAnnealingModel{T} 
    c0p::T = -63.37 - 0.212log(3600) # Yamada et al. 2007 zircon, converted from log(h) to log(s)
    c1p::T = 0.212      # Yamada et al. 2007 zircon
    bp::T = 43.00       # Yamada et al. 2007 zircon
    l0::T = 11.17       # [um] effective initial track length (μmax)
    l0_sigma::T = 0.051 # [um] effective initial track length uncertainty (σ)
end

"""
```julia
Guenthner2013FC(
    C0 = 6.24534         # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    C1 = -0.11977        # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    C2 = -314.93688      # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    C3 = -14.2868        # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    alpha = -0.057206897 # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    l0 = 11.17           # [um] Initial track length
    l0_sigma = 0.051     # [um] Initial track length uncertainty
)
```
Fanning Curvilinear zircon annealing model with simplified Box-Cox transform 
from Guenthner et al. 2013 (doi: 10.2475/03.2013.01)
"""
Base.@kwdef struct Guenthner2013FC{T<:AbstractFloat} <: FanningCurvilinearZircon{T} 
    C0::T = 6.24534         # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    C1::T = -0.11977        # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    C2::T = -314.93688      # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    C3::T = -14.2868        # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    alpha::T = -0.057206897 # Guenthner et al. 2013 re-fit of Yamada et al. 2007 zircon
    l0::T = 11.17           # [um] Initial track length
    l0_sigma::T = 0.051     # [um] Initial track length uncertainty
end

"""
```julia
Ketcham1999FC(
    C0::T = -19.844     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    C1::T = 0.38951     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    C2::T = -51.253     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    C3::T = -7.6423     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    alpha::T = -0.12327 # Box-Cox transform parameter
    beta::T = -11.988   # Box-Cox transform parameter
    l0::T = 16.38       # [um] Initial track length
    l0_sigma::T = 0.09  # [um] Initial track length unertainty
)
```
Fanning Curvilinear apatite annealing model from Ketcham, 1999 (doi: 10.2138/am-1999-0903)
"""
Base.@kwdef struct Ketcham1999FC{T<:AbstractFloat} <: FanningCurvilinearApatite{T} 
    C0::T = -19.844     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    C1::T = 0.38951     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    C2::T = -51.253     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    C3::T = -7.6423     # "Simultaneous fit" from Ketcham et al. 1999 apatite
    alpha::T = -0.12327 # Box-Cox transform parameter
    beta::T = -11.988   # Box-Cox transform parameter
    l0::T = 16.38       # [um] Initial track length
    l0_sigma::T = 0.09  # [um] Initial track length unertainty
end

"""
```julia
Ketcham2007FC(
    C0 = 0.39528     # "Simultaneous fit" from Ketcham et al. 2007 apatite
    C1 = 0.01073     # "Simultaneous fit" from Ketcham et al. 2007 apatite
    C2 = -65.12969   # "Simultaneous fit" from Ketcham et al. 2007 apatite
    C3 = -7.91715    # "Simultaneous fit" from Ketcham et al. 2007 apatite
    alpha = 0.04672  # "Simultaneous fit" from Ketcham et al. 2007 apatite
    l0 = 16.38       # [um] Initial track length
    l0_sigma = 0.09  # [um] Initial track length unertainty
)
```
Fanning Curvilinear apatite annealing model with simplified Box-Cox transform 
from Ketcham, 2007 (doi: 10.2138/am.2007.2281)
"""
Base.@kwdef struct Ketcham2007FC{T<:AbstractFloat} <: FanningCurvilinearApatite{T} 
    C0::T = 0.39528     # "Simultaneous fit" from Ketcham et al. 2007 apatite
    C1::T = 0.01073     # "Simultaneous fit" from Ketcham et al. 2007 apatite
    C2::T = -65.12969   # "Simultaneous fit" from Ketcham et al. 2007 apatite
    C3::T = -7.91715    # "Simultaneous fit" from Ketcham et al. 2007 apatite
    alpha::T = 0.04672  # "Simultaneous fit" from Ketcham et al. 2007 apatite
    l0::T = 16.38       # [um] Initial track length
    l0_sigma::T = 0.09  # [um] Initial track length unertainty
end

"""
```julia
Jones2021FA(
    C0 = 1.374           # Annealing parameter
    C1 = -4.192812e-5    # Annealing parameter
    C2 = -22.70885029    # Annealing parameter
    C3 = 0.0             # Annealing parameter
    l0 = 10.60           # [um] Initial track length
    l0_sigma = 0.19      # [um] Initial track length uncertainty
)
```
Parallel Arrhenius monazite annealing model modified from
Jones et al., 2021 (doi: 10.5194/gchron-3-89-2021)
"""
Base.@kwdef struct Jones2021FA{T<:AbstractFloat} <: MonaziteAnnealingModel{T} 
    C0::T = 1.374           # Annealing parameter
    C1::T = -4.192812e-5    # Annealing parameter
    C2::T = -22.70885029    # Annealing parameter
    C3::T = 0.0             # Annealing parameter
    l0::T = 10.60           # [um] Initial track length
    l0_sigma::T = 0.19      # [um] Initial track length uncertainty
end

## --- Define DiffusivityModel types

abstract type DiffusivityModel{T} end
abstract type ZirconHeliumModel{T} <: DiffusivityModel{T} end
abstract type ApatiteHeliumModel{T} <: DiffusivityModel{T} end

# Implement methods to allow broadcasting
Base.length(x::DiffusivityModel) = 1
Base.iterate(x::DiffusivityModel) = (x, nothing)
Base.iterate(x::DiffusivityModel, state) = nothing

"""
```julia
ZRDAAM(
    DzD0::T = 193188.0          # Diffusivity [cm^2/sec], crystalline endmember
    DzD0_logsigma::T=1/2        # log units (default = 1/2 = a factor of ℯ two-sigma)
    DzEa::T=165.0               # Activation energy [kJ/mol], crystalline endmember
    DzEa_logsigma::T=1/2        # log units (default = 1/2 = a factor of ℯ two-sigma)
    DN17D0::T = 6.367E-3        # Diffusivity [cm^2/sec], amorphous endmember
    DN17D0_logsigma::T=1/2      # log units (default = 1/2 = a factor of ℯ two-sigma)
    DN17Ea::T=71.0              # Activation energy [kJ/mol], amorphous endmember
    DN17Ea_logsigma::T=1/2      # log units (default = 1/2 = a factor of ℯ two-sigma)
    lint0::T=45920.0            # [nm]
    SV::T=1.669                 # [1/nm]
    Bα::T=5.48E-19              # Amorphous material produced per alpha decay [g/alpha]
    Phi::T=3.0                  # unitless
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
    DzD0::T = 193188.0          # Diffusivity [cm^2/sec], crystalline endmember
    DzD0_logsigma::T=1/2        # log units (default = 1/2 = a factor of ℯ two-sigma)
    DzEa::T=165.0               # Activation energy [kJ/mol], crystalline endmember
    DzEa_logsigma::T=1/2        # log units (default = 1/2 = a factor of ℯ two-sigma)
    DN17D0::T = 6.367E-3        # Diffusivity [cm^2/sec], amorphous endmember
    DN17D0_logsigma::T=1/2      # log units (default = 1/2 = a factor of ℯ two-sigma)
    DN17Ea::T=71.0              # Activation energy [kJ/mol], amorphous endmember
    DN17Ea_logsigma::T=1/2      # log units (default = 1/2 = a factor of ℯ two-sigma)
    lint0::T=45920.0            # [nm]
    SV::T=1.669                 # [1/nm]
    Bα::T=5.48E-19              # Amorphous material produced per alpha decay [g/alpha]
    Phi::T=3.0                  # unitless
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
    D0L::T=0.6071               # Diffusivity [cm^2/s]
    D0L_logsigma::T=1/2         # log units (default = 1/2 = a factor of ℯ two-sigma)
    EaL::T=122.3                # Activation energy [kJ/mol]
    EaL_logsigma::T=1/2         # log units (default = 1/2 = a factor of ℯ two-sigma)
    EaTrap::T=34.0              # Activation energy [kJ/mol]
    EaTrap_logsigma::T=1/2      # log units (default = 1/2 = a factor of ℯ two-sigma)
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
    D0L::T=0.6071               # Diffusivity [cm^2/s]
    D0L_logsigma::T=1/2         # log units (default = 1/2 = a factor of ℯ two-sigma)
    EaL::T=122.3                # Activation energy [kJ/mol]
    EaL_logsigma::T=1/2         # log units (default = 1/2 = a factor of ℯ two-sigma)
    EaTrap::T=34.0              # Activation energy [kJ/mol]
    EaTrap_logsigma::T=1/2      # log units (default = 1/2 = a factor of ℯ two-sigma)
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
    rmr0_sigma::T=0.15          # Damage conversion parameter uncertainty
    kappa::T=1.04-0.83          # Damage conversion parameter
    kappa_rmr0::T=1.04          # Damage conversion parameter (the sum of kappa and rmr0)
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
    admdist::Vector{<:ApatiteHeliumModel{T}}
    zdmdist::Vector{<:ZirconHeliumModel{T}}
end
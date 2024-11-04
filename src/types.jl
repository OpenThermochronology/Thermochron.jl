## --- Define DamageModel types
abstract type DamageModel end
export DamageModel

Base.@kwdef struct ZRDAAM{T<:AbstractFloat} <: DamageModel 
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
    rmr0::T=0.2                 # Zircon anealing parameter (not normally called this in ZRDAAM context, but effectively the same thing)
end
export ZRDAAM

Base.@kwdef struct RDAAM{T<:AbstractFloat} <: DamageModel 
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
    kappa::T=1.04-0.83          # Damage conversion parameter
    kappa_rmr0::T=1.04          # Damage conversion parameter (the sum of kappa and rmr0)
end
export RDAAM

## --- Define Boundary type to specify the working area
struct Boundary{T<:AbstractFloat}
    agepoints::Vector{T} # Ma
    Tpoints::Vector{T}   # Degrees C
    Tpointsₚ::Vector{T}   # Degrees C
    T₀::Vector{T}
    ΔT::Vector{T}
    npoints::Int
end
function Boundary(T::Type=Float64; agepoints, Tpoints, T₀, ΔT)
    @assert length(agepoints) == length(Tpoints) == length(T₀) == length(ΔT)
    Boundary{T}(T.(agepoints), 
        T.(Tpoints), 
        T.(Tpoints), 
        T.(T₀), 
        T.(ΔT), 
        length(agepoints)
    )
end
export Boundary

## -- Define Constraint type to specify aribitrary t-T constraint boxes
struct Constraint{T<:AbstractFloat}
    agepoints::Vector{T} # Ma
    Tpoints::Vector{T}   # Degrees C
    agepointsₚ::Vector{T} # Ma
    Tpointsₚ::Vector{T}   # Degrees C
    Age₀::Vector{T}
    ΔAge::Vector{T}
    T₀::Vector{T}
    ΔT::Vector{T}
    npoints::Int
end
function Constraint(T::Type=Float64; agepoints=Float64[], Tpoints=Float64[], Age₀=Float64[], ΔAge=Float64[], T₀=Float64[], ΔT=Float64[])
    agepoints = isempty(agepoints) ? Age₀ + ΔAge/2 : agepoints
    Tpoints = isempty(Tpoints) ? T₀ + ΔT/2 : Tpoints
    @assert length(agepoints) == length(Tpoints) == length(Age₀) == length(ΔAge) == length(T₀) == length(ΔT)
    Constraint{T}(T.(agepoints),
        T.(Tpoints),
        T.(agepoints), 
        T.(Tpoints), 
        T.(Age₀), 
        T.(ΔAge), 
        T.(T₀), 
        T.(ΔT), 
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


struct Boundary{T<:AbstractFloat}
    agePoints::Vector{T} # Ma
    TPoints::Vector{T}   # Degrees C
    TPointsₚ::Vector{T}   # Degrees C
    T₀::Vector{T}
    ΔT::Vector{T}
    npoints::Int
end
function Boundary(T::Type=Float64; agePoints, TPoints, T₀, ΔT)
    agePoints = T.(agePoints)
    TPoints =  T.(TPoints)
    T₀ = T.(T₀)
    ΔT = T.(ΔT)
    @assert length(agePoints) == length(TPoints) == length(T₀) == length(ΔT)
    Boundary(agePoints, TPoints, copy(TPoints), T₀, ΔT, length(agePoints))
end
export Boundary

struct Unconformity{T<:AbstractFloat}
    agePoints::Vector{T} # Ma
    TPoints::Vector{T}   # Degrees C
    agePointsₚ::Vector{T} # Ma
    TPointsₚ::Vector{T}   # Degrees C
    Age₀::Vector{T}
    ΔAge::Vector{T}
    T₀::Vector{T}
    ΔT::Vector{T}
    npoints::Int
end
function Unconformity(T::Type=Float64; agePoints=Float64[], TPoints=Float64[], Age₀=Float64[], ΔAge=Float64[], T₀=Float64[], ΔT=Float64[])
    agePoints = T.(agePoints)
    TPoints = T.(TPoints)
    Age₀ = T.(Age₀)
    ΔAge = T.(ΔAge)
    T₀ = T.(T₀)
    ΔT = T.(ΔT)
    @assert length(agePoints) == length(TPoints) == length(Age₀) == length(ΔAge) == length(T₀) == length(ΔT)
    Unconformity(agePoints, TPoints, copy(agePoints), copy(TPoints), Age₀, ΔAge, T₀, ΔT, length(agePoints))
end
export Unconformity

struct DetailInterval{T<:AbstractFloat}
    agemin::T
    agemax::T
    minpoints::Int
end
function DetailInterval(T::Type=Float64; agemin, agemax, minpoints)
    DetailInterval(T(agemin), T(agemax), Int(minpoints))
end
export DetailInterval

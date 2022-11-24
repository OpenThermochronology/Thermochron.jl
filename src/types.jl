
struct Boundary{T<:AbstractFloat}
    agepoints::Vector{T} # Ma
    Tpoints::Vector{T}   # Degrees C
    Tpointsₚ::Vector{T}   # Degrees C
    T₀::Vector{T}
    ΔT::Vector{T}
    npoints::Int
end
function Boundary(T::Type=Float64; agepoints, Tpoints, T₀, ΔT)
    agepoints = T.(agepoints)
    Tpoints =  T.(Tpoints)
    T₀ = T.(T₀)
    ΔT = T.(ΔT)
    @assert length(agepoints) == length(Tpoints) == length(T₀) == length(ΔT)
    Boundary(agepoints, Tpoints, copy(Tpoints), T₀, ΔT, length(agepoints))
end
export Boundary

struct Unconformity{T<:AbstractFloat}
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
function Unconformity(T::Type=Float64; agepoints=Float64[], Tpoints=Float64[], Age₀=Float64[], ΔAge=Float64[], T₀=Float64[], ΔT=Float64[])
    agepoints = T.(agepoints)
    Tpoints = T.(Tpoints)
    Age₀ = T.(Age₀)
    ΔAge = T.(ΔAge)
    T₀ = T.(T₀)
    ΔT = T.(ΔT)
    @assert length(agepoints) == length(Tpoints) == length(Age₀) == length(ΔAge) == length(T₀) == length(ΔT)
    Unconformity(agepoints, Tpoints, copy(agepoints), copy(Tpoints), Age₀, ΔAge, T₀, ΔT, length(agepoints))
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

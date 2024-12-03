# Compact show for chronometers
function Base.show(io::IO, x::T) where {T<:HeliumSample}
    t = Base.typename(T).wrapper
    σ = round(x.age_sigma, sigdigits=2)
    d = floor(Int, log10(x.age)-log10(x.age_sigma))
    μ = round(x.age, sigdigits=2+d)
    print(io, "$t($(μ)±$(σ) Ma)")
end
function Base.show(io::IO, x::T) where {T<:FissionTrackSample}
    t = Base.typename(T).wrapper
    σ = round(x.age_sigma, sigdigits=2)
    d = floor(Int, log10(x.age)-log10(x.age_sigma))
    μ = round(x.age, sigdigits=2+d)
    print(io, "$t($(μ)±$(σ) Ma)")
end
function Base.show(io::IO, x::T) where {T<:FissionTrackLength}
    t = Base.typename(T).wrapper
    l = round(x.length, sigdigits=3)
    print(io, "$t($(l) μm)")
end

# Compact show for annealing models
function Base.show(io::IO, x::T) where {T<:AnnealingModel}
    if get(io, :compact, false)
        print(io, "$T(…)")
    else
        Base.show_default(io, x)
    end
end

# Compact show for diffusivity models
function Base.show(io::IO, x::T) where {T<:DiffusivityModel}
    if get(io, :compact, false)
        print(io, "$T(…)")
    else
        Base.show_default(io, x)
    end
end
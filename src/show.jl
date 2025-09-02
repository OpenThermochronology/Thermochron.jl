# Compact show for chronometers
function Base.show(io::IO, x::T) where {T<:MultipleDomain}
    t = Base.typename(T).wrapper
    n = length(x.domains)
    c = Base.typename(typeof(first(x.domains))).wrapper
    print(io, "$t($n $c)")
end
function Base.show(io::IO, x::T) where {T<:ArgonSample}
    t = Base.typename(T).wrapper
    σ = round(x.age_sigma, sigdigits=2)
    d = log10(x.age)-log10(x.age_sigma)
    μ = round(x.age, sigdigits=2+floor(Int, d*!isnan(d)))
    print(io, "$t($(μ)±$(σ) Ma)")
end
function Base.show(io::IO, x::T) where {T<:HeliumSample}
    t = Base.typename(T).wrapper
    σ = round(x.age_sigma, sigdigits=2)
    d = log10(x.age)-log10(x.age_sigma)
    μ = round(x.age, sigdigits=2+floor(Int, d*!isnan(d)))
    print(io, "$t($(μ)±$(σ) Ma)")
end
function Base.show(io::IO, x::T) where {T<:FissionTrackSample}
    t = Base.typename(T).wrapper
    σ = round(x.age_sigma, sigdigits=2)
    d = log10(x.age)-log10(x.age_sigma)
    μ = round(x.age, sigdigits=2+floor(Int, d*!isnan(d)))
    print(io, "$t($(μ)±$(σ) Ma)")
end
function Base.show(io::IO, x::T) where {T<:FissionTrackLength}
    t = Base.typename(T).wrapper
    l = round(x.length, sigdigits=3)
    print(io, "$t($(l) μm)")
end
function Base.show(io::IO, x::T) where {T<:ApatiteTrackLengthOriented}
    t = Base.typename(T).wrapper
    l = round(x.length, sigdigits=3)
    θ = round(x.angle, sigdigits=3)
    print(io, "$t($(l) μm, $(θ)°)")
end

# Verbose show methods
printshort(x::AbstractArray) = "[$(first(x)) … $(last(x))]"
printshort(x::AbstractRange) = "$x"
function Base.show(io::IO, ::MIME"text/plain", x::T) where {T<:MultipleDomain}
    print(io, """$T:
      step age          : $(printshort(x.step_age)) Ma
      step age sigma    : $(printshort(x.step_age_sigma)) Ma
      fraction degassed : $(printshort(x.fraction_experimental))
      offset            : $(x.offset) C from the surface
      domains           : $(printshort(x.domains))
      volume fraction   : $(printshort(x.volume_fraction))
      agesteps          : $(printshort(agediscretization(x))) Ma
    """
    )
end
function Base.show(io::IO, ::MIME"text/plain", x::T) where {T<:ArgonSample}
    print(io, """$T:
      age       : $(x.age) Ma
      age_sigma : $(x.age_sigma) Ma
      offset    : $(x.offset) C from the surface
      K-40      : $(printshort(x.r40K / (6.022E23 / 1E6 / 39.96399848))) ppm
      rsteps    : $(x.rsteps) μm
      agesteps  : $(x.agesteps) Ma
    """
    )
end
function Base.show(io::IO, ::MIME"text/plain", x::T) where {T<:HeliumSample}
    print(io, """$T:
      age       : $(x.age) Ma
      age_sigma : $(x.age_sigma) Ma
      offset    : $(x.offset) C from the surface
      U-238     : $(printshort(x.r238U / (6.022E23 / 1E6 / 238))) ppm
      U-235     : $(printshort(x.r235U / (6.022E23 / 1E6 / 235))) ppm
      Th-232    : $(printshort(x.r232Th / (6.022E23 / 1E6 / 232))) ppm
      Sm-147    : $(printshort(x.r147Sm / (6.022E23 / 1E6 / 147))) ppm
      rsteps    : $(printshort(x.rsteps)) μm
      agesteps  : $(printshort(x.agesteps)) Ma
    """
    )
end
function Base.show(io::IO, ::MIME"text/plain", x::T) where {T<:FissionTrackSample}
    print(io, """$T:
      age       : $(x.age) Ma
      age_sigma : $(x.age_sigma) Ma
      offset    : $(x.offset) C from the surface
      agesteps  : $(x.agesteps) Ma
    """
    )
end
function Base.show(io::IO, ::MIME"text/plain", x::T) where {T<:ApatiteFT}
    print(io, """$T:
      age       : $(x.age) Ma
      age_sigma : $(x.age_sigma) Ma
      offset    : $(x.offset) C from the surface
      rmr0      : $(x.rmr0)
      agesteps  : $(printshort(x.agesteps)) Ma
    """
    )
end
function Base.show(io::IO, ::MIME"text/plain", x::T) where {T<:FissionTrackLength}
    print(io, """$T:
      length    : $(x.length) μm
      offset    : $(x.offset) C from the surface
      agesteps  : $(printshort(x.agesteps)) Ma
    """
    )
end
function Base.show(io::IO, ::MIME"text/plain", x::T) where {T<:ApatiteTrackLengthOriented}
    print(io, """$T:
      length    : $(x.length) μm
      angle     : $(x.angle) degrees from c-axis
      offset    : $(x.offset) C from the surface
      rmr0      : $(x.rmr0)
      agesteps  : $(printshort(x.agesteps)) Ma
    """
    )
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
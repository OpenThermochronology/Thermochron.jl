## ---  Various useful internal utility functions

    # This is not public API and so may be liable to change in future Julia versions
    # but can seemingly mean a factor of two in performance over using `step` function 
    # according to benchmarks based on use of `step_at` in this package
    @inline step_fast(x::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}) = x.step.hi + x.step.lo 
    @inline step_fast(x::StepRangeLen) = step(x) # All other StepRangeLen
    @inline step_fast(x::OrdinalRange) = step(x)

    # Get local step size for an AbstractVector at a given index
    function step_at(x::AbstractVector{T}, i::Integer) where {T<:Number}
        length(x) > 1 || return zero(T)
        f,l = firstindex(x), lastindex(x)
        if f < i < l
            # Assume boundaries are halfway between nodes, implying a step width 
            # of (x[i+1]-x[i])/2 + (x[i]-x[i-1])/2, or simplified
            T((x[i+1] - x[i-1])/2)
        elseif i <= f
            # Assume first step is symmetric
            (x[f+1] - x[f])
        else
            # Assume last step is symmetric
            (x[l] - x[l-1])
        end
    end
    @inline step_at(x::Union{OrdinalRange, StepRangeLen}, ::Integer) = step_fast(x)

    function leftbound_at(x::AbstractVector{T}, i::Integer) where {T<:Number}
        Tf = float(T)
        length(x) > 1 || Tf(NaN)
        f,l = firstindex(x), lastindex(x)
        # Assume boundaries are halfway between nodes
        if f < i < l
            Tf((x[i-1] + x[i])/2)
        elseif i <= f
            Tf(x[f] - (x[f+1] - x[f])/2)
        else
            Tf((x[l-1] + x[l])/2)
        end
    end
    function rightbound_at(x::AbstractVector{T}, i::Integer) where {T<:Number}
        Tf = float(T)
        length(x) > 1 || Tf(NaN)
        f,l = firstindex(x), lastindex(x)
        # Assume boundaries are halfway between nodes
        if f < i < l
            Tf((x[i] + x[i+1])/2)
        elseif i <= f
            Tf((x[f] + x[f+1])/2)
        else
            Tf(x[l] + (x[l] - x[l-1])/2)
        end
    end


    # Normal log likelihood
    @inline function norm_ll(mu::Number, sigma::Number, x::Number)
        δ = x - mu
        σ² = sigma^2
        -0.5*(log(2*pi*σ²) + δ^2/σ²)
    end
    @inline function norm_ll(mu::Number, sigma::Number, x::Number, x_sigma::Number)
        δ = x - mu
        σ² = sigma^2 + x_sigma^2
        -0.5*(log(2*pi*σ²) + δ^2/σ²)
    end

    function copyunique!(dest, source)
        isempty(source) && return dest
        id = firstindex(dest)
        for i in eachindex(source)
            sᵢ = source[i]
            add = true
            for i in firstindex(dest):id
                if dest[i] == sᵢ
                    add=false
                    break
                end
            end
            if add
                dest[id] = sᵢ
                id += 1
            end
        end
        return dest
    end

    function volumefraction(shape::Symbol, redges::AbstractArray, r::Number=last(redges))
        if shape === :spherical
            return (redges[2:end].^3 .- redges[1:end-1].^3)./r^3
        elseif shape === :cylindrical
            return (redges[2:end].^2 .- redges[1:end-1].^2)./r^2
        elseif shape === :planar
            return (redges[2:end] .- redges[1:end-1])./r
        else
            @error "Shape $shape not recognized; options are `:spherical`, `:cylindrical`, or `:planar`"
        end
    end

    """
    ```julia
    sphereintersectionfraction(r₁, r₂, d)
    ```
    Calculate the fraction of the surface area of a sphere s₂ with radius `r₂`
    that intersects the interior of a sphere s₁ of radius `r₁` if the two are
    separated by distance `d`.
    """
    function sphereintersectionfraction(r₁::T, r₂::T, d::T) where T <: AbstractFloat
        # Let r₁ and r₂ be the radii of two spheres s₁ and s₂, separated by distance d
        @assert (r₁ >= 0) &&  (r₂ >= 0)
        d = abs(d)
        if d > r₁+r₂ ## If separated by more than dmax, there is no intersection
            omega = zero(T)
        elseif d > abs(r₁-r₂)
            # X is the radial distance between the center of s2 and the interseciton plane
            # See e.g., http://mathworld.wolfram.com/Sphere-SphereIntersection.html
            # x = (d^2 - r₁^2 + r₂^2) / (2 * d)

            # Let omega be the solid angle of intersection normalized by 4pi,
            # where the solid angle of a cone is 2pi*(1-cos(theta)) and cos(theta)
            # is adjacent/hypotenuse = x/r₂.
            # omega = (1 - x/r₂)/2

            # Rearranged for optimization:
            omega = one(T)/2 - (d^2 - r₁^2 + r₂^2) / (4 * d * r₂)

        elseif r₁<r₂ # If r₁ is entirely within r₂
            omega = zero(T)
        else
            omega = one(T) # If r₂ is entirely within r₁
        end
        return omega
    end

    """
    ```julia
    slabsphereintersectionfraction(rₚ, rₛ, d)
    ```
    Calculate the fraction of the surface area of a sphere s with radius `rₛ`
    that intersects the interior of a planar slab p of halfwidth `rₚ` if the two are
    separated by distance `d`.
    """
    function slabsphereintersectionfraction(rₚ::T, rₛ::T, d::T) where T <: AbstractFloat
        # Let rₚ be the halfwidth of a planar slab and rₛ the radius of a sphere
        # separated by distance d
        @assert (rₚ >= 0) && (rₛ >= 0)
        d = abs(d)
        if d > (rₚ+rₛ) ## If separated by more than dmax, there is no intersection
            omega = zero(T)
        elseif d > abs(rₚ-rₛ)
            # X is the radial distance between the center of the sphere s and the interseciton plane
            x = d - rₚ

            # Let omega be the solid angle of intersection normalized by 4pi,
            # where the solid angle of a cone is 2pi*(1-cos(theta)) and cos(theta)
            # is adjacent/hypotenuse = x/rₛ.
            omega = (one(T) - x/rₛ)/2  # zero when x = rₛ, one when x = -rₛ

        elseif rₚ < rₛ 
            # If halfwidth of the planar slab is less than radius of the sphere
            # then both sides of the slab are intersecting the sphere 
            
            # # Larger solid angle of intersection
            # x₊ = d - rₚ
            # omega₊ = (one(T) - x₊/rₛ)/2 
            # # Smaller solid angle of intersection 
            # x₋ = d + rₚ
            # omega₋ = (one(T) - x₋/rₛ)/2  
            # # Difference = intersection within slab
            # omega = omega₊ - omega₋

            # Rearranged and terms cancelled
            omega = rₚ/rₛ
        else
            # If sphere is entirely within planar slab
            omega = one(T) 
        end
        return omega
    end


    """
    ```julia
    sphereintersectiondensity(redges::Vector, rvolumes::Vector, ralpha, d)
    ```
    Calculate the volume-nomalized fractional intersection density of an alpha
    stopping sphere of radius `ralpha` with each concentric shell (with shell edges
    `redges` and relative volumes `rvolumes`) of a spherical crystal where the
    two are separated by distance `d`
    """
    function sphereintersectiondensity(redges::AbstractVector{T}, rvolumes::AbstractVector{T}, ralpha::T, d::T) where T <: AbstractFloat
        sphereintersectiondensity!(zeros(T, length(redges)-1), redges, rvolumes, ralpha, d)
    end
    function sphereintersectiondensity!(dint::DenseVector{T}, redges::AbstractVector{T}, rvolumes::AbstractVector{T}, ralpha::T, d::T) where T <: AbstractFloat
        @assert eachindex(dint) == eachindex(rvolumes) == firstindex(redges):lastindex(redges)-1
        fintlast = sphereintersectionfraction(first(redges), ralpha, d)
        @inbounds for i ∈ eachindex(dint,rvolumes)
            # Integrated intersection fraction for each concentric sphere (redges) of crystal
            fint = sphereintersectionfraction(redges[i+1], ralpha, d)

            # Intersection fraction for each spherical shell of the crystal (subtracting
            # one concentric sphere from the next) normalized by shell volume
            dint[i] = (fint-fintlast) / rvolumes[i]
            fintlast = fint
        end
        return dint
    end

    """
    ```julia
    slabsphereintersectiondensity(redges::Vector, ralpha, d)
    ```
    Calculate the fractional intersection density of an alpha
    stopping sphere of radius `ralpha` with each concentric slab-shell (with edges
    `redges` and relative volumes `rvolumes`) of a planar slab crystal where the
    two are separated by distance `d`
    """
    function slabsphereintersectiondensity(redges::AbstractVector{T}, ralpha::T, d::T) where T <: AbstractFloat
        slabsphereintersectiondensity!(zeros(T, length(redges)-1), redges, ralpha, d)
    end
    function slabsphereintersectiondensity!(dint::DenseVector{T}, redges::AbstractVector{T}, ralpha::T, d::T) where T <: AbstractFloat
        @assert eachindex(dint) == firstindex(redges):lastindex(redges)-1
        fintlast = slabsphereintersectionfraction(first(redges), ralpha, d)
        @inbounds for i ∈ eachindex(dint)
            # Integrated intersection fraction for each concentric sphere (redges) of crystal
            fint = slabsphereintersectionfraction(redges[i+1], ralpha, d)

            # Intersection fraction for each slab-shell of the crystal (subtracting
            # one concentric slab from the next)
            dint[i] = fint-fintlast
            fintlast = fint
        end
        return dint
    end

    """
    ```julia
    alphastoppingpower(mineral)
    ```
    Alpha particle stopping power relative to that of apatite,
    as calculated from the mean alpha stopping distances of  
    Ketcham et al. 2011 (doi: 10.1016/j.gca.2011.10.011)
    """
    alphastoppingpower(mineral::String) = alphastoppingpower(Symbol(lowercase(mineral)))
    function alphastoppingpower(mineral::Symbol)
        if mineral===:apatite
            1.0
        elseif mineral===:zircon
            1.2082235609872582
        elseif mineral===:titanite || mineral===:sphene
            1.076593548244311
        elseif mineral===:monazite   
            1.1633813225022829
        elseif mineral===:xenotime   
            1.2369421530297124
        elseif mineral===:rutile     
            1.227818975881152
        elseif mineral===:magnetite  
            1.3482551879811657
        elseif mineral===:hematite   
            1.3860091775480607
        elseif mineral===:goethite   
            1.2106968510177099
        elseif mineral===:barite     
            1.0358153480322465
        else
            @warn "Mineral $mineral not recognized.
            Alpha stopping power unknown; using average."
            1.189 # Average of all the above
        end
    end
    export alphastoppingpower


    """
    ```julia
    alphacorrectionspherical(alpharadii, rparent, redges, relvolumes, parentₒ, redgesₒ, relvolumesₒ)
    ```
    Effective parent concentrations accounting for alphas per decay at secular
    equilibrium as well as alpha ejection (given internal radial `rparent` 
    concentration vector and spatial discretization `redges`, `relvolumes`) and
    injection (given uniform external `parentₒ` concentration and spatial 
    discretization `redgesₒ`, `relvolumesₒ`) for an infinite slab mineral grain.
    """
    function alphacorrectionspherical(alpharadii, rparent::AbstractVector{T}, redges, relvolumes, parentₒ::T, redgesₒ, relvolumesₒ) where {T}
        # Prepare
        rsteps = cntr(redges)
        rstepsₒ = cntr(redgesₒ)
        dint = zeros(T, length(rsteps))
        parentHe = zeros(T, length(rparent))
        @assert eachindex(rparent) == eachindex(rsteps) == eachindex(relvolumes) 
        @assert eachindex(rparent) == eachindex(dint) == eachindex(parentHe)
        @assert eachindex(rstepsₒ) == eachindex(relvolumesₒ)

        # Correct for alpha ejection
        @inbounds for i in eachindex(rsteps, relvolumes, rparent)
            for rₐ in alpharadii
                sphereintersectiondensity!(dint,redges,relvolumes,rₐ,rsteps[i])
                @. parentHe += relvolumes[i] * dint * rparent[i]
            end
        end
        # Correct for alpha injection
        if parentₒ > 0
            @inbounds for i in eachindex(rstepsₒ, relvolumesₒ)
                for rₐ in alpharadii
                    (rstepsₒ[i] - first(rstepsₒ)) > rₐ && continue
                    sphereintersectiondensity!(dint, redges, relvolumes, rₐ, rstepsₒ[i])
                    @. parentHe += relvolumesₒ[i] * dint * parentₒ
                end
            end
        end
        return parentHe
    end

    """
    ```julia
    alphacorrectionslab(alpharadii, rparent, redges, parentₒ, redgesₒ)
    ```
    Effective parent concentrations accounting for alphas per decay at secular
    equilibrium as well as alpha ejection (given internal radial `rparent` 
    concentration vector and spatial discretization `redges`) and injection
    (given uniform external `parentₒ` concentration and spatial discretization
    `redgesₒ`) for an infinite slab mineral grain.
    """
    function alphacorrectionslab(alpharadii, rparent::AbstractVector{T}, redges, parentₒ::T, redgesₒ) where {T}
        # Prepare
        rsteps = cntr(redges)
        rstepsₒ = cntr(redgesₒ)
        dint = zeros(T, length(rsteps))
        parentHe = zeros(T, length(rparent))
        @assert eachindex(rparent) == eachindex(rsteps)
        @assert eachindex(rparent) == eachindex(dint) == eachindex(parentHe)

        # Correct for alpha ejection
        @inbounds for i in eachindex(rsteps, rparent)
            for alpharadius in alpharadii
                slabsphereintersectiondensity!(dint, redges,alpharadius, rsteps[i])
                @. parentHe += dint * rparent[i]
            end
        end
        # Correct for alpha injection
        if parentₒ > 0
            @inbounds for i in eachindex(rstepsₒ)
                for alpharadius in alpharadii
                    (rstepsₒ[i] - first(rstepsₒ)) > alpharadius && continue
                    slabsphereintersectiondensity!(dint, redges, alpharadius, rstepsₒ[i])
                    @. parentHe += dint * parentₒ
                end
            end
        end
        return parentHe
    end


    # Utility function for agepoint and Tpoint buffers
    function collectto!(buffer, a, b, c)
        i₀ = firstindex(buffer)
        copyto!(buffer, i₀, a, 1, length(a))
        i₀ += length(a)
        copyto!(buffer, i₀, b, 1, length(b))
        i₀ += length(b)
        copyto!(buffer, i₀, c, 1, length(c))
        i₀ += length(c)
        return view(buffer, firstindex(buffer):i₀-1)
    end

    # Utility functions for checking the number of distinct t-T nodes in a given time interval
    function pointsininterval(points::AbstractArray, npoints::Int, lower::Number, upper::Number)
        @assert firstindex(points) == 1
        @assert npoints <= lastindex(points)
        n = 0
        @inbounds for i = 1:npoints
            if  lower < points[i] <= upper
                n += 1
            end
        end
        return n
    end
    function pointsininterval(points::AbstractArray, npoints::Int, lower::Number, upper::Number, nodes::AbstractArray)
        @assert firstindex(points) == 1
        @assert npoints <= lastindex(points)
        n = 0
        @inbounds for j in eachindex(nodes)
            if lower < nodes[j] <= upper
                for k in 1:npoints
                    l, u = minmax(leftbound_at(nodes, j), rightbound_at(nodes, j))
                    if max(lower,l) <= points[k] < min(upper, u)
                        # Only one point per interval can count towards total
                        n += 1
                        break
                    end
                end
            end
        end
        return n
    end
    proposeddetail(path::TtPath, npointsₚ::Int, detail::DetailInterval) = pointsininterval(path.agepointsₚ, npointsₚ, detail.agemin, detail.agemax, path.agesteps)
    accepteddetail(path::TtPath, npoints::Int, detail::DetailInterval) = pointsininterval(path.agepoints, npoints, detail.agemin, detail.agemax, path.agesteps)
    mindetail(path::TtPath, npoints::Int, detail::DetailInterval) = min(accepteddetail(path, npoints, detail), detail.minpoints)

    # Check if point k is distinct from other points in list within ± δ
    function isdistinct(points::AbstractArray, k::Int, δ::Number, npoints::Int=length(points))
        @assert k ∈ eachindex(points)
        notused = max(length(points)-npoints, 0)
        @inbounds for i in Iterators.drop(reverse(eachindex(points)), notused)
            if i!=k && abs(points[i] - points[k]) <= δ
                return false
            end
        end
        return true
    end
    function isdistinct(points::AbstractArray, k::Int, nodes::AbstractArray, npoints::Int=length(points); kwargs...)
        length(points) > 1 || return false
        i = searchsortedfirst(nodes, points[k]; kwargs...)
        i₋,i₊ = firstindex(nodes), lastindex(nodes)
        if i₋ < i <= i₊
            l, u = nodes[i-1], nodes[i]
        elseif i <= i₋
            l, u = typemin(nodes[i₋]), nodes[i₋]
        else
            l, u = nodes[i₊], typemax(nodes[i₊])
        end
        @assert k ∈ eachindex(points)
        notused = max(length(points)-npoints, 0)
        @inbounds for i in Iterators.drop(reverse(eachindex(points)), notused)
            if (i != k) && (l < points[i] <= u)
                return false
            end
        end
        return true
    end


    # Utility functions for checking maximum reheating or cooling rate
    function maxdiff(x::AbstractArray{T}) where {T}
        i₀ = firstindex(x)
        δₘ = zero(T)
        if length(x) > 1
            last = x[i₀]
            @inbounds for i ∈ (i₀+1):(i₀+length(x)-1)
                δᵢ = x[i] - last
                if δᵢ > δₘ
                    δₘ = δᵢ
                end
                last = x[i]
            end
        end
        return δₘ
    end
    function mindiff(x::AbstractArray{T}) where {T}
        i₀ = firstindex(x)
        δₘ = zero(T)
        if length(x) > 1
            last = x[i₀]
            @inbounds for i ∈ (i₀+1):(i₀+length(x)-1)
                δᵢ = x[i] - last
                if δᵢ < δₘ
                    δₘ = δᵢ
                end
                last = x[i]
            end
        end
        return δₘ
    end
    function maxabsdiff(x::AbstractArray{T}) where {T}
        i₀ = firstindex(x)
        δₘ = zero(T) 
        if length(x) > 1
            last = x[i₀]
            @inbounds for i ∈ (i₀+1):(i₀+length(x)-1)
                δᵢ = abs(x[i] - last)
                if δᵢ > δₘ
                    δₘ = δᵢ
                end
                last = x[i]
            end
        end
        return δₘ
    end

    function diff_ll(x::AbstractArray{T}, μ::Number, σ::Number; maxsigma::Number=6) where {T<:Number}
        i₀ = firstindex(x)
        d = Normal(μ, σ)
        ll = zero(float(T))
        if length(x) > 1
            last = x[i₀]
            @inbounds for i ∈ (i₀+1):lastindex(x)
                δᵢ = x[i] - last
                if δᵢ > 0
                    ll += logccdf(d, δᵢ)
                    if δᵢ > (μ + maxsigma*σ)
                        # Firm cap at some number standard deviations above the mean
                        ll -= maxintfloat(T)
                    end
                end
                last = x[i]
            end
        end
        return ll
    end

    # Log likelihood functions for ZRDAAM and RDAAM kinetic uncertainties
        
    function kinetic_ll(dmₚ::ZRDAAM, dm::ZRDAAM)
        norm_ll(log(dm.DzD0), dm.DzD0_logsigma, log(dmₚ.DzD0)) + 
        norm_ll(log(dm.DzEa), dm.DzEa_logsigma, log(dmₚ.DzEa)) + 
        norm_ll(log(dm.DN17D0), dm.DN17D0_logsigma, log(dmₚ.DN17D0)) + 
        norm_ll(log(dm.DN17Ea), dm.DN17Ea_logsigma, log(dmₚ.DN17Ea)) +
        norm_ll(dm.rmin, dm.rmin_sigma, dmₚ.rmin)
    end
    function kinetic_ll(dmₚ::RDAAM, dm::RDAAM)
        norm_ll(log(dm.D0L), dm.D0L_logsigma, log(dmₚ.D0L)) + 
        norm_ll(log(dm.EaL), dm.EaL_logsigma, log(dmₚ.EaL)) + 
        norm_ll(log(dm.EaTrap), dm.EaTrap_logsigma, log(dmₚ.EaTrap))+
        norm_ll(dm.rmr0, dm.rmr0_sigma, dmₚ.rmr0)
    end
    function kinetic_ll(dmₚ::Diffusivity, dm::Diffusivity)
        norm_ll(log(dm.D0), dm.D0_logsigma, log(dmₚ.D0)) + 
        norm_ll(log(dm.Ea), dm.Ea_logsigma, log(dmₚ.Ea))
    end
    function kinetic_ll(dmₚ::MDDiffusivity{T}, dm::MDDiffusivity{T}) where {T<:AbstractFloat}
        ll = zero(T)
        for i in eachindex(dm.Ea, dm.Ea_logsigma, dm.D0, dm.D0_logsigma)
            ll += norm_ll(log(dm.Ea[i]), dm.Ea_logsigma[i], log(dmₚ.Ea[i]))
            ll += norm_ll(log(dm.D0[i]), dm.D0_logsigma[i], log(dmₚ.D0[i]))
        end
        return ll
    end
    function kinetic_ll(dmₚ::SDDiffusivity, dm::SDDiffusivity)
        kinetic_ll(unwrap(dmₚ), unwrap(dm)) + 
        norm_ll(log(dm.scale), dm.scale_logsigma, log(dmₚ.scale))
    end
    function kinetic_ll(dmₚ::AnnealingModel{T}, dm::AnnealingModel{T}) where {T<:AbstractFloat}
        @assert dmₚ == dm
        return zero(T)
    end
    function kinetic_ll!(updatekinetics::BitVector, damodelsₚ::Vector{<:Model{T}}, damodels::Vector{<:Model{T}}) where {T}
        fill!(updatekinetics, true)
        ll = zero(T)
        # Count each unique kinetic model exactly once towards the kinetic log likelihood
        for i in eachindex(damodelsₚ, damodels)
            if updatekinetics[i]
                dm = damodels[i]
                dmₚ = damodelsₚ[i]
                ll += kinetic_ll(dmₚ, dm)
                for j in eachindex(damodelsₚ, damodels)
                    if updatekinetics[j] && (damodels[j] == dm) && (damodelsₚ[i] == dmₚ)
                        updatekinetics[j] = false
                    end
                end
            end
        end
        return ll
    end


    function draw_from_population(track::FissionTrackLength{T}, bandwidth::T) where {T<:AbstractFloat}
        pr = track.pr::Vector{T}
        rΣ = sum(pr)*rand()
        rΣ > 0 || return T(NaN)
        Σ = zero(T)
        i = firstindex(pr)
        while i < lastindex(pr)
            Σ += pr[i]
            if Σ < rΣ
                i += 1
            else
                break
            end
        end
        return rand(Normal{T}(track.r[i], bandwidth))
    end
    function draw_from_population(x::AbstractVector, cumulative_fraction::AbstractVector)
        @assert eachindex(x) == eachindex(cumulative_fraction)
        i = 0
        r = last(cumulative_fraction) * rand()
        for f in cumulative_fraction
            i += 1
            f > r && break
        end
        return x[i]
    end
    
    """
    ```julia
    simannealT(n::Integer, Tₐ::Number, λₐ::Number)
    ```
    To avoid getting stuck in local optima, increase the probability 
    of accepting new proposals at higher annealing "temperature"

    """
    function simannealT(n::Integer, Tₐ::Number, λₐ::Number)
        @assert λₐ >= 0
        @assert Tₐ >= 0
        return Tₐ * exp(-λₐ*n) + 1
    end


    # Utitlity functions for dealing with boundary conditions
    function reflecting(x, xmin, xmax)
        @assert xmin < xmax
        Δx = xmax-xmin
        d,r = divrem(x-xmin, Δx)
        if isodd(d)
            xmax - abs(r)
        else
            xmin + abs(r)
        end
    end
    function periodic(x, xmin, xmax)
        @assert xmin < xmax
        Δx = xmax-xmin
        xmin + mod(x-xmin, Δx)
    end
    function hard(x, xmin, xmax)
        @assert xmin < xmax
        min(max(x, xmin), xmax)
    end

    textrema(path::TtPath) = textrema(path.boundary)
    function textrema(boundary::Boundary)
        a, b = first(boundary.agepoints), last(boundary.agepoints)
        a < b ? (a, b) : (b, a)
    end
    Textrema(path::TtPath) = Textrema(path.boundary)
    function Textrema(boundary::Boundary)
        a, b = first(boundary.T₀), last(boundary.T₀)
        a < b ? (a, b) : (b, a)
    end

    function boundtime(t::Number, boundary::Boundary)
        tmin, tmax = textrema(boundary)
        if boundary.tboundary === :reflecting
            reflecting(t, tmin, tmax)
        elseif boundary.tboundary === :hard
            hard(t, tmin, tmax)
        elseif boundary.tboundary === :periodic
            periodic(t, tmin, tmax)
        else
            @warn "`tboundary` $(boundary.tboundary) not recognized; choose `:reflecting`, `:hard`, or `:periodic`.\nDefaulting to `:reflecting`."
            reflecting(t, tmin, tmax)
        end
    end
    function boundtemp(T::Number, boundary::Boundary)
        Tmin, Tmax = Textrema(boundary)
        if boundary.Tboundary === :reflecting
            reflecting(T, Tmin, Tmax)
        elseif boundary.Tboundary === :hard
            hard(T, Tmin, Tmax)
        elseif boundary.Tboundary === :periodic
            periodic(T, Tmin, Tmax)
        else
            @warn "`Tboundary` $(boundary.Tboundary) not recognized; choose `:reflecting`, `:hard`, or `:periodic`.\nDefaulting to `:reflecting`."
            reflecting(T, Tmin, Tmax)
        end
    end

    # Move a t-T point and apply boundary conditions
    function movepoint!(path::TtPath{T}, k::Int, npoints::Int) where {T}
        # Move the age of one model point
        path.agepointsₚ[k] += rand(Normal{T}(zero(T), path.σⱼtₚ[k]))

        # Move the Temperature of one model point
        path.Tpointsₚ[k] += rand(Normal{T}(zero(T), path.σⱼTₚ[k]))

        # Apply time boundary conditions
        path.agepointsₚ[k] = boundtime(path.agepointsₚ[k], path.boundary)

        # Apply Temperature boundary conditions
        path.Tpointsₚ[k] = boundtemp(path.Tpointsₚ[k], path.boundary)

        # Ensure uniqueness
        while !isdistinct(path.agepointsₚ, k, path.agesteps, npoints; rev=true)
            path.agepointsₚ[k] += rand(Normal{T}(zero(T), path.σⱼtₚ[k]))
            path.agepointsₚ[k] = boundtime(path.agepointsₚ[k], path.boundary)
        end

        return path
    end

    # Add a t-T point
    function addpoint!(path::TtPath{T}, k::Int, npoints::Int) where {T}
        @assert eachindex(path.agepointsₚ) == eachindex(path.Tpointsₚ) == eachindex(path.σⱼtₚ) == eachindex(path.σⱼTₚ)

        tmin, tmax = textrema(path.boundary)
        Tmin, Tmax = Textrema(path.boundary)

        # Pick an age uniformly from the model age domain
        path.agepointsₚ[k] = rand(path.agesteps)

        # Find the closest existing points (if any)
        ages = view(path.agepointsₚ, Base.OneTo(k-1))
        i₋ = findclosestbelow(path.agepointsₚ[k], ages)
        i₊ = findclosestabove(path.agepointsₚ[k], ages)

        # Find values for the closest younger point
        inbounds₋ = firstindex(ages) <= i₋ <= lastindex(ages)
        t₋ = inbounds₋ ? path.agepointsₚ[i₋] : tmin
        T₋ = inbounds₋ ? path.Tpointsₚ[i₋] : Tmin
        σⱼt₋ = inbounds₋ ? path.σⱼtₚ[i₋] : (tmax-tmin)/60
        σⱼT₋ = inbounds₋ ? path.σⱼTₚ[i₋] : (Tmax-Tmin)/60

        # Find values for the closest older point
        inbounds₊ = firstindex(ages) <= i₊ <= lastindex(ages)
        t₊ = inbounds₊ ? path.agepointsₚ[i₊] : tmax
        T₊ = inbounds₊ ? path.Tpointsₚ[i₊] : Tmax
        σⱼt₊ = inbounds₊ ? path.σⱼtₚ[i₊] : (tmax-tmin)/60
        σⱼT₊ = inbounds₊ ? path.σⱼTₚ[i₊] : (Tmax-Tmin)/60

        # Interpolate
        f = (path.agepointsₚ[k] - t₋) / (t₊ - t₋)
        f *= !isnan(f)
        path.Tpointsₚ[k] = f*T₊ + (1-f)*T₋
        path.σⱼtₚ[k] = f*σⱼt₊ + (1-f)*σⱼt₋
        path.σⱼTₚ[k] = f*σⱼT₊ + (1-f)*σⱼT₋

        # Move the point from the interpolated value
        return movepoint!(path, k, npoints)
    end

    function replacepoint!(path::TtPath{T}, k::Int, n::Int) where {T}
        path.agepointsₚ[k] = path.agepointsₚ[n]
        path.Tpointsₚ[k] = path.Tpointsₚ[n]
        path.σⱼtₚ[k] = path.σⱼtₚ[n]
        path.σⱼTₚ[k] = path.σⱼTₚ[n]
        return path
    end

    # Adjust initial and final t-T boundaries
    function movebounds!(boundary::Boundary)
        @inbounds for i in eachindex(boundary.Tpointsₚ)
            boundary.Tpointsₚ[i] = boundary.T₀[i] + rand()*boundary.ΔT[i]
        end
        return boundary
    end
    function movebounds!(constraint::Constraint, boundary::Boundary)
        @inbounds for i in eachindex(constraint.agepointsₚ, constraint.Tpointsₚ)
            constraint.agepointsₚ[i] = boundtime(rand(constraint.agedist[i]), boundary)
            constraint.Tpointsₚ[i] = boundtemp(rand(constraint.Tdist[i]), boundary)
        end
        return constraint
    end
    function movebounds!(path::TtPath)
        movebounds!(path.constraint, path.boundary)
        movebounds!(path.boundary)
        return path
    end
    
    function randomize!(boundary::Boundary)
        @inbounds for i in eachindex(boundary.Tpoints)
            boundary.Tpoints[i] = boundary.T₀[i] + rand()*boundary.ΔT[i]
        end
        return boundary
    end
    function randomize!(constraint::Constraint, boundary::Boundary)
        @inbounds for i in eachindex(constraint.agepoints, constraint.Tpoints)
            constraint.agepoints[i] = boundtime(rand(constraint.agedist[i]), boundary)
            constraint.Tpoints[i] = boundtemp(rand(constraint.Tdist[i]), boundary)
        end
        return constraint
    end

    function initialproposal!(path::TtPath, npoints::Int)
        randomize!(path.boundary)
        randomize!(path.constraint, path.boundary)
        collectaccepted!(path, 0)
        for i in eachindex(path.agepoints, path.Tpoints)
            if (i-firstindex(path.agepoints)) < path.detail.minpoints
                path.agepoints[i] = rand(Uniform(path.detail.agemin, path.detail.agemax)) # Uniformly from the detail interval
            else
                path.agepoints[i] = rand(path.agesteps) # Uniformly from the model age domain
            end
            path.Tpoints[i] = linterp1(reverse(path.agesteps), path.Tsteps, first(path.agesteps)-path.agepoints[i])
        end
        collectaccepted!(path, npoints)
        return path
    end

    function collectaccepted!(path::TtPath{T}, npoints::Int) where {T<:AbstractFloat}
        agepoints = view(path.agepoints, Base.OneTo(npoints))
        ages = collectto!(path.agepointbuffer, agepoints, path.boundary.agepoints, path.constraint.agepoints)
        Tpoints = view(path.Tpoints, Base.OneTo(npoints))
        temperatures = collectto!(path.Tpointbuffer, Tpoints, path.boundary.Tpoints, path.constraint.Tpoints)
        linterp1s!(path.Tsteps, path.knot_index, ages, temperatures, path.agesteps)
        return path
    end
    function collectproposal!(path::TtPath{T}, npointsₚ::Int) where {T<:AbstractFloat}
        agepointsₚ = view(path.agepointsₚ, Base.OneTo(npointsₚ))
        ages = collectto!(path.agepointbuffer, agepointsₚ, path.boundary.agepoints, path.constraint.agepointsₚ)
        Tpointsₚ = view(path.Tpointsₚ, Base.OneTo(npointsₚ))
        temperatures = collectto!(path.Tpointbuffer, Tpointsₚ, path.boundary.Tpointsₚ, path.constraint.Tpointsₚ)
        linterp1s!(path.Tsteps, path.knot_index, ages, temperatures, path.agesteps)
        return path
    end
    function resetproposal!(path::TtPath)
        copyto!(path.agepointsₚ, path.agepoints)
        copyto!(path.Tpointsₚ, path.Tpoints)
        copyto!(path.σⱼtₚ, path.σⱼt)
        copyto!(path.σⱼTₚ, path.σⱼT)
        copyto!(path.constraint.agepointsₚ, path.constraint.agepoints)
        copyto!(path.constraint.Tpointsₚ, path.constraint.Tpoints)
        copyto!(path.boundary.Tpointsₚ, path.boundary.Tpoints)
    end
    function acceptproposal!(path::TtPath)
        copyto!(path.agepoints, path.agepointsₚ)
        copyto!(path.Tpoints, path.Tpointsₚ)
        copyto!(path.σⱼt, path.σⱼtₚ)
        copyto!(path.σⱼT, path.σⱼTₚ)
        copyto!(path.constraint.agepoints, path.constraint.agepointsₚ)
        copyto!(path.constraint.Tpoints, path.constraint.Tpointsₚ)
        copyto!(path.boundary.Tpoints, path.boundary.Tpointsₚ)
    end

    # Update jumping distribution
    function updatejumping!(path::TtPath, k::Int)
        if path.agepointsₚ[k] != path.agepoints[k]
            path.σⱼtₚ[k] = ℯ * abs(path.agepointsₚ[k] - path.agepoints[k])
        end
        if path.Tpointsₚ[k] != path.Tpoints[k]
            path.σⱼTₚ[k] = ℯ * abs(path.Tpointsₚ[k] - path.Tpoints[k])
        end
        return path
    end

    # Utilities for dealing with kinetic models that wrap other kinetic models
    unwrap(x::Model) = x
    unwrap(x::SDDiffusivity) = x.model
    wrap(x::Model, ::Model) = x
    wrap(x::DiffusivityModel, y::SDDiffusivity) = SDDiffusivity(x, y.scale, y.scale_logsigma)

    # Adjust kinetic models
    movekinetics(dm::Model, p=0.5) = dm
    function movekinetics(zdm::ZRDAAM{T}, p=0.5) where {T}
        ZRDAAM(
            DzEa = (rand()<p) ? exp(log(zdm.DzEa)+randn(T)*zdm.DzEa_logsigma/2) : zdm.DzEa,
            DzEa_logsigma = zdm.DzEa_logsigma,
            DzD0 = (rand()<p) ? exp(log(zdm.DzD0)+randn(T)*zdm.DzD0_logsigma/2) : zdm.DzD0,
            DzD0_logsigma = zdm.DzD0_logsigma,
            DN17Ea = (rand()<p) ? exp(log(zdm.DN17Ea)+randn(T)*zdm.DN17Ea_logsigma/2) : zdm.DN17Ea,
            DN17Ea_logsigma = zdm.DN17Ea_logsigma,
            DN17D0 = (rand()<p) ? exp(log(zdm.DN17D0)+randn(T)*zdm.DN17D0_logsigma/2) : zdm.DN17D0,
            DN17D0_logsigma = zdm.DN17D0_logsigma,
            rmin = (rand()<p) ? reflecting(zdm.rmin + randn(T)*zdm.rmin_sigma/2, 0, 1) : zdm.rmin,
            rmin_sigma = zdm.rmin_sigma,
        )
    end
    function movekinetics(adm::RDAAM{T}, p=0.5) where {T}
        rmr0 = (rand()<p) ? reflecting(adm.rmr0 + randn(T)*adm.rmr0_sigma/2, 0, 1) : adm.rmr0
        RDAAM(
            D0L = (rand()<p) ? exp(log(adm.D0L)+randn(T)*adm.D0L_logsigma/2) : adm.D0L,
            D0L_logsigma = adm.D0L_logsigma,
            EaL = (rand()<p) ? exp(log(adm.EaL)+randn(T)*adm.EaL_logsigma/2) : adm.EaL,
            EaL_logsigma = adm.EaL_logsigma,
            EaTrap = (rand()<p) ? exp(log(adm.EaTrap)+randn(T)*adm.EaTrap_logsigma/2) : adm.EaTrap,
            EaTrap_logsigma = adm.EaTrap_logsigma,
            rmr0 = rmr0,
            rmr0_sigma = adm.rmr0_sigma,
            kappa = adm.kappa_rmr0 - rmr0,
        )
    end
    function movekinetics(dm::Diffusivity{T}, p=0.5) where {T}
        Diffusivity(
            D0 = (rand()<p) ? exp(log(dm.D0)+randn(T)*dm.D0_logsigma/2) : dm.D0,
            D0_logsigma = dm.D0_logsigma,
            Ea = (rand()<p) ? exp(log(dm.Ea)+randn(T)*dm.Ea_logsigma/2) : dm.Ea,
            Ea_logsigma = dm.Ea_logsigma,
        )
    end
    function movekinetics(dm::MDDiffusivity{T}, p=0.5) where {T}
        MDDiffusivity(
            D0 = ((rand()<p) ? @.(exp(log(dm.D0)+(rand()<p)*randn(T)*dm.D0_logsigma/4)) : dm.D0),
            D0_logsigma = dm.D0_logsigma,
            Ea = ((rand()<p) ? @.(exp(log(dm.Ea)+(rand()<p)*randn(T)*dm.Ea_logsigma/10)) : dm.Ea),
            Ea_logsigma = dm.Ea_logsigma,
        )
    end
    function movekinetics(dm::SDDiffusivity{T}, p=0.5) where {T}
        SDDiffusivity(
            model = movekinetics(dm.model, p),
            scale = (rand()<p) ? exp(log(dm.scale)+randn(T)*dm.scale_logsigma/2) : dm.scale,
            scale_logsigma = dm.scale_logsigma,
        )
    end
    movewrapperkinetics(dm::Model, p=0.5) = dm
    function movewrapperkinetics(dm::SDDiffusivity{T}, p=0.5) where {T}
        SDDiffusivity(
            model = dm.model,
            scale = (rand()<p) ? exp(log(dm.scale)+randn(T)*dm.scale_logsigma/2) : dm.scale,
            scale_logsigma = dm.scale_logsigma,
        )
    end
    function movekinetics!(damodels::Vector{<:Model}, updatekinetics::BitVector, p=0.5)
        fill!(updatekinetics, true)
        for i in eachindex(damodels)
            if updatekinetics[i]
                dm = damodels[i]
                dmₚ = movekinetics(dm, p)
                # Keep diffusivity models which start identical, identical
                for j in eachindex(damodels)
                    if updatekinetics[j]
                        if damodels[j] == dm
                            # Identical models
                            damodels[j] = dmₚ
                            updatekinetics[j] = false
                        elseif damodels[j] == unwrap(dm)
                            # Model that we wrap
                            damodels[j] = unwrap(dmₚ)
                            updatekinetics[j] = false 
                        elseif unwrap(damodels[j]) == unwrap(dm)
                            # Model that wraps us (or that wraps same thing as us)
                            damodels[j] = wrap(unwrap(dmₚ), movewrapperkinetics(damodels[j], p))
                            updatekinetics[j] = false 
                        end
                    end
                end
            end
        end
        return damodels
    end

    # Adjust model uncertainties of chronometers
    function movesigma!(σcalc::AbstractVector{T}, chrons::AbstractVector{<:Chronometer}) where {T<:AbstractFloat}
        for C in (ZirconFT, MonaziteFT, ApatiteFT, ZirconHe, ApatiteHe, SphericalHe, PlanarHe, SphericalAr, PlanarAr, MultipleDomain)
            r = abs(randn(T))
            for i in eachindex(σcalc, chrons)
                if chrons[i] isa C
                    σcalc[i] *= r
                end
            end
        end
        return σcalc
    end

    # Utility functions to calculate model ages for all chronometers at once
    function model(chrons::Vector{<:Chronometer{T}}, damodels::Vector{<:Model{T}}, Tsteps::AbstractVector{T}; kwargs...) where {T}
        μcalc = zeros(T, size(chrons))
        σcalc =  zeros(T, size(chrons))
        ll = model!(μcalc, σcalc, chrons, damodels, Tsteps; kwargs...)
        return μcalc, σcalc, ll
    end
    function model!(μcalc::AbstractVector{T}, σcalc::AbstractVector{T}, chrons::Vector{<:Chronometer{T}}, damodels::Vector{<:Model{T}}, Tsteps::AbstractVector{T}; 
            rescale::Bool=false,
            rescalestepheating::Bool=true,  
            redegastracer::Bool=false, 
            stepwisetracerfraction::Bool=false,
            partitiondaughter::Bool=false, 
        ) where {T<:AbstractFloat}
        imax = argmax(i->length(timediscretization(chrons[i])), eachindex(chrons))
        tsteps = timediscretization(chrons[imax])
        @assert issorted(tsteps)
        if eachindex(tsteps) != eachindex(Tsteps)
            @info chrons[imax]
            @info tsteps
            @info Tsteps
        end
        @assert eachindex(tsteps) == eachindex(Tsteps) "`tsteps` has indices $(eachindex(tsteps)) while `Tsteps` has indices $(eachindex(Tsteps))"

        # Pre-anneal ZRDAAM samples, if any
        if any(x-> eltype(x) <: ZRDAAM, damodels)
            zdm = unwrap(damodels[findfirst(x-> eltype(x) <: ZRDAAM, damodels)])::ZRDAAM{T}
            anneal!(chrons, ZirconHe, tsteps, Tsteps, zdm)
        end
        # Pre-anneal RDAAM samples, if any
        if any(x-> eltype(x) <: RDAAM, damodels)
            adm = unwrap(damodels[findfirst(x-> eltype(x) <: RDAAM, damodels)])::RDAAM{T}
            anneal!(chrons, ApatiteHe, tsteps, Tsteps, adm)
        end

        # Flag any fission track lengths for recalculation
        for i in eachindex(chrons)
            chrons[i] isa FissionTrackLength && (chrons[i].calc .= NaN)
        end

        # Optionally rescale log likelihoods to avoid one chronometer type from dominating the inversion
        scalegar = rescale ? sqrt(count(x->(isa(x, SphericalAr)||isa(x, PlanarAr)), chrons)) : 1
        scaleghe = rescale ? sqrt(count(x->(isa(x, SphericalHe)||isa(x, PlanarHe)), chrons)) : 1
        scalezhe = rescale ? sqrt(count(x->isa(x, ZirconHe), chrons)) : 1
        scaleahe = rescale ? sqrt(count(x->isa(x, ApatiteHe), chrons)) : 1
        scalezft = rescale ? sqrt(count(x->isa(x, ZirconFT), chrons)) : 1
        scalemft = rescale ? sqrt(count(x->isa(x, MonaziteFT), chrons)) : 1
        scaleaft = rescale ? sqrt(count(x->isa(x, ApatiteFT), chrons)) : 1
        scaleztl = rescale ? sqrt(count(x->isa(x, ZirconTrackLength), chrons)) : 1
        scalemtl = rescale ? sqrt(count(x->isa(x, MonaziteTrackLength), chrons)) : 1
        scaleatl = rescale ? sqrt(count(x->isa(x, ApatiteTrackLength), chrons)) : 1
        scaleato = rescale ? sqrt(count(x->isa(x, ApatiteTrackLengthOriented), chrons)) : 1
        scalesdd = rescale ? sqrt(count(x->isa(x, SingleDomain), chrons)) : 1
        scalemdd = rescale ? sqrt(count(x->isa(x, MultipleDomain), chrons)) : 1

        # Cycle through each Chronometer, model and calculate log likelihood
        ll = zero(T)
        for i in eachindex(chrons, μcalc, σcalc, damodels)
            c, dm = chrons[i], damodels[i]
            first_index = firstindex(Tsteps) + length(tsteps) - length(timediscretization(c))
            Tstepsᵢ = @views(Tsteps[first_index:end])
            if isa(c, SphericalAr) || isa(c, PlanarAr)
                c::Union{SphericalAr{T}, PlanarAr{T}}
                μcalc[i] = modelage(c, Tstepsᵢ, dm::Diffusivity{T}; partitiondaughter)
                ll += norm_ll(μcalc[i], σcalc[i], value(c), stdev(c))/scalegar
            elseif isa(c, SphericalHe) || isa(c, PlanarHe)
                c::Union{SphericalHe{T}, PlanarHe{T}}
                μcalc[i] = modelage(c, Tstepsᵢ, dm::Diffusivity{T}; partitiondaughter)
                ll += norm_ll(μcalc[i], σcalc[i], value(c), stdev(c))/scaleghe
            elseif isa(c, ZirconHe)
                c::ZirconHe{T}
                μcalc[i] = modelage(c, Tstepsᵢ, dm::ZirconHeliumModel{T}; partitiondaughter)
                ll += norm_ll(μcalc[i], σcalc[i], value(c), stdev(c))/scalezhe
            elseif isa(c, ApatiteHe)
                c::ApatiteHe{T}
                μcalc[i] = modelage(c, Tstepsᵢ, dm::ApatiteHeliumModel{T}; partitiondaughter)
                ll += norm_ll(μcalc[i], σcalc[i], value(c), stdev(c))/scaleahe
            elseif isa(c, ZirconFT)
                c::ZirconFT{T}
                μcalc[i] = modelage(c, Tstepsᵢ, dm::ZirconAnnealingModel{T})
                ll += norm_ll(μcalc[i], σcalc[i], value(c), stdev(c))/scalezft
            elseif isa(c, MonaziteFT)
                c::MonaziteFT{T}
                μcalc[i] = modelage(c, Tstepsᵢ, dm::MonaziteAnnealingModel{T})
                ll += norm_ll(μcalc[i], σcalc[i], value(c), stdev(c))/scalemft
            elseif isa(c, ApatiteFT)
                c::ApatiteFT{T}
                μcalc[i] = modelage(c, Tstepsᵢ, dm::ApatiteAnnealingModel{T})
                ll += norm_ll(μcalc[i], σcalc[i], value(c), stdev(c))/scaleaft
            elseif isa(c, ZirconTrackLength)
                c::ZirconTrackLength{T}
                if isnan(first(c.calc))
                    μ, σcalc[i] = modellength(c, Tstepsᵢ, dm::ZirconAnnealingModel{T})
                else
                    μ, σcalc[i] = c.calc
                end
                μcalc[i] = max(draw_from_population(c, σcalc[i]), zero(T)) # Draw from dist, but cannot be negative
                ll += model_ll(c, σcalc[i])/scaleztl
            elseif isa(c, MonaziteTrackLength)
                c::MonaziteTrackLength{T}
                if isnan(first(c.calc))
                    μ, σcalc[i] = modellength(c, Tstepsᵢ, dm::MonaziteAnnealingModel{T})
                else
                    μ, σcalc[i] = c.calc
                end
                μcalc[i] = max(draw_from_population(c, σcalc[i]), zero(T)) # Draw from dist, but cannot be negative
                ll += model_ll(c, σcalc[i])/scalemtl
            elseif isa(c, ApatiteTrackLength)
                c::ApatiteTrackLength{T}
                if isnan(first(c.calc))
                    μ, σcalc[i] = modellength(c, Tstepsᵢ, dm::ApatiteAnnealingModel{T})
                else
                    μ, σcalc[i] = c.calc
                end
                μcalc[i] = max(draw_from_population(c, σcalc[i]), zero(T)) # Draw from dist, but cannot be negative
                ll += model_ll(c, σcalc[i])/scaleatl
            elseif isa(c, ApatiteTrackLengthOriented)
                c::ApatiteTrackLengthOriented{T}
                if isnan(first(c.calc))
                    μ, σcalc[i] = modellength(c, Tstepsᵢ, dm::ApatiteAnnealingModel{T})
                else
                    μ, σcalc[i] = c.calc
                end
                μcalc[i] = max(draw_from_population(c, σcalc[i]), zero(T)) # Draw from dist, but cannot be negative
                ll += model_ll(c, σcalc[i])/scaleato
            elseif isa(c, SingleDomain)
                c::SingleDomain{T, <:Union{ArgonSample{T}, HeliumSample{T}}}
                age, stepage, fraction = modelage(c, Tstepsᵢ, dm::DiffusivityModel{T}; redegastracer, partitiondaughter)
                @assert issorted(fraction) "Degassing fraction is not properly cumulative"
                if redegastracer
                    if stepwisetracerfraction
                        ll += stepwise_degassing_ll(c; rescale=rescalestepheating)/scalesdd
                    else
                        ll += cumulative_degassing_ll(c; rescale=rescalestepheating)/scalesdd
                    end
                end
                μcalc[i] = draw_from_population(stepage, fraction)
                ll += model_ll(c, σcalc[i]; rescale=rescalestepheating)/scalesdd
            elseif isa(c, MultipleDomain)
                c::MultipleDomain{T, <:Union{PlanarAr{T}, SphericalAr{T}}}
                stepage, fraction = modelage(c, Tstepsᵢ, dm::MDDiffusivity{T}; redegastracer, partitiondaughter)
                @assert issorted(fraction) "Degassing fraction is not properly cumulative"
                if redegastracer
                    if stepwisetracerfraction
                        ll += stepwise_degassing_ll(c; rescale=rescalestepheating)/scalemdd
                    else
                        ll += cumulative_degassing_ll(c; rescale=rescalestepheating)/scalemdd
                    end
                end
                μcalc[i] = draw_from_population(stepage, fraction)
                ll += model_ll(c, σcalc[i]; rescale=rescalestepheating)/scalemdd
            else
                # NaN if not calculated
                μcalc[i] = T(NaN)
            end
        end
            
        return ll
    end

## --- Extracting things from KineticResults

# Query Diffusivities and MDDiffusivities from a KineticResult
function Diffusivity(kr::KineticResult)
    any(x->isa(x, Diffusivity), kr) || return nothing
    ia = findall(x->isa(x, Diffusivity), kr[:,1])
    return collect(kr[ia,:]')
end
function MDDiffusivity(kr::KineticResult)
    any(x->isa(x, MDDiffusivity), kr) || return nothing
    ia = findall(x->isa(x, MDDiffusivity), kr[:,1])
    return collect(kr[ia,:]')
end

## --- End of File
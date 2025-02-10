## ---  Various useful internal utility functions

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

    """
    ```julia
    intersectionfraction(r₁, r₂, d)
    ```
    Calculate the fraction of the surface area of a sphere s2 with radius `r₂`
    that intersects the interior of a sphere s1 of radius `r₁` if the two are
    separated by distance `d`. Assumptions: `r₁`>0, `r₂`>0, `d`>=0
    """
    function intersectionfraction(r₁::T, r₂::T, d::T) where T <: AbstractFloat
        dmax = r₁+r₂
        dmin = abs(r₁-r₂)
        if d > dmax ## If separated by more than dmax, there is no intersection
            omega = zero(T)
        elseif d > dmin
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
    intersectiondensity(redges::Vector, rvolumes::Vector, ralpha, d)
    ```
    Calculate the volume-nomalized fractional intersection density of an alpha
    stopping sphere of radius `ralpha` with each concentric shell (with shell edges
    `redges` and relative volumes `rvolumes`) of a spherical crystal where the
    two are separated by distance `d`
    """
    function intersectiondensity(redges::AbstractVector{T}, rvolumes::AbstractVector{T}, ralpha::T, d::T) where T <: AbstractFloat
        dint = Array{T}(undef, length(redges) - 1)
        intersectiondensity!(dint, redges, rvolumes, ralpha, d)
    end
    function intersectiondensity!(dint::DenseVector{T}, redges::AbstractVector{T}, rvolumes::AbstractVector{T}, ralpha::T, d::T) where T <: AbstractFloat
        I = firstindex(redges):lastindex(redges)-1
        @assert eachindex(dint) == eachindex(rvolumes) == I
        fintlast = intersectionfraction(first(redges), ralpha, d)
        @inbounds for i ∈ I
            # Integrated intersection fraction for each concentric sphere (redges) of crystal
            fint = intersectionfraction(redges[i+1], ralpha, d)

            # Intersection fraction for each spherical shell of the crystal (subtracting
            # one concentric sphere from the next) normalized by shell volume
            dint[i] = (fint-fintlast) / rvolumes[i]
            fintlast = fint
        end
        return dint
    end

    # Alpha particle stopping power relative to apatite
    # as calculated from the mean alpha stopping distances of  
    # Ketcham et al. 2011 (doi: 10.1016/j.gca.2011.10.011)
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

    # Utility functions for checking the nummber of distinct t-T nodes in a given time interval
    function pointsininterval(points::DenseArray, npoints::Int, min::Number, max::Number)
        n = 0
        @inbounds for i = 1:npoints
            if  min < points[i] < max
                n += 1
            end
        end
        return n
    end
    function pointsininterval(points::DenseArray, npoints::Int, min::Number, max::Number, δ::Number)
        n = 0
        @inbounds for i = 1:npoints
            if  min < points[i] < max
                n += isdistinct(points, i, δ, npoints)
            end
        end
        return n
    end

    # Check if point k is distinct from other points in list within ± δ
    function isdistinct(points::AbstractArray, k::Int, δ::Number, npoints::Int=length(points))
        @assert npoints <= length(points)
        I = firstindex(points):firstindex(points)+npoints-1
        @inbounds for i in I
            if i!=k && abs(points[i] - points[k]) <= δ
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

    function diff_ll(x::AbstractArray, μ::Number, σ::Number)
        i₀ = firstindex(x)
        inv_s2 = 1/(2*σ*σ)
        ll = zero(typeof(inv_s2))
        if length(x) > 1
            last = x[i₀]
            @inbounds for i ∈ (i₀+1):(i₀+length(x)-1)
                δᵢ = x[i] - last
                if δᵢ > μ
                    δμ = δᵢ - μ
                    ll -= δμ * δμ * inv_s2
                end
                last = x[i]
            end
        end
        return ll
    end

    # Log likelihood functions for ZRDAAM and RDAAM kinetic uncertainties
    function kinetic_ll(zdmₚ::ZRDAAM, zdm::ZRDAAM)
        norm_ll(log(zdm.DzD0), zdm.DzD0_logsigma, log(zdmₚ.DzD0)) + 
        norm_ll(log(zdm.DzEa), zdm.DzEa_logsigma, log(zdmₚ.DzEa)) + 
        norm_ll(log(zdm.DN17D0), zdm.DN17D0_logsigma, log(zdmₚ.DN17D0)) + 
        norm_ll(log(zdm.DN17Ea), zdm.DN17Ea_logsigma, log(zdmₚ.DN17Ea)) +
        norm_ll(zdm.rmr0, zdm.rmr0_sigma, zdmₚ.rmr0)
    end
    function kinetic_ll(admₚ::RDAAM, adm::RDAAM)
        norm_ll(log(adm.D0L), adm.D0L_logsigma, log(admₚ.D0L)) + 
        norm_ll(log(adm.EaL), adm.EaL_logsigma, log(admₚ.EaL)) + 
        norm_ll(log(adm.EaTrap), adm.EaTrap_logsigma, log(admₚ.EaTrap))+
        norm_ll(adm.rmr0, adm.rmr0_sigma, admₚ.rmr0)
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
    function movepoint!(path::TtPath{T}, k::Int) where {T}
        # Move the age of one model point
        path.agepointsₚ[k] += rand(Normal{T}(zero(T), path.σⱼtₚ[k]))

        # Move the Temperature of one model point
        path.Tpointsₚ[k] += rand(Normal{T}(zero(T), path.σⱼTₚ[k]))

        # Apply time boundary conditions
        path.agepointsₚ[k] = boundtime(path.agepointsₚ[k], path.boundary)

        # Apply Temperature boundary conditions
        path.Tpointsₚ[k] = boundtemp(path.Tpointsₚ[k], path.boundary)

        return path
    end

    # Add a t-T point
    function addpoint!(path::TtPath{T}, k::Int) where {T}
        @assert eachindex(path.agepointsₚ) == eachindex(path.Tpointsₚ) == eachindex(path.σⱼtₚ) == eachindex(path.σⱼTₚ)

        tmin, tmax = textrema(path.boundary)
        Tmin, Tmax = Textrema(path.boundary)

        # Pick an age uniformly within the boundaries
        path.agepointsₚ[k] = rand(Uniform(tmin, tmax))

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
        return movepoint!(path, k)
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
    function randomize!(agepoints::AbstractVector, Tpoints::AbstractVector, boundary::Boundary, detail::DetailInterval)
        agemin, agemax = textrema(boundary)
        Tmin, Tmax = Textrema(boundary)
        for i in eachindex(agepoints, Tpoints)
            if (i-firstindex(agepoints)) < detail.minpoints
                agepoints[i] = rand(Uniform(detail.agemin, detail.agemax))
            else
                agepoints[i] = rand(Uniform(agemin, agemax))
            end
            Tpoints[i] = rand(Uniform(Tmin, Tmax))
        end
        return agepoints, Tpoints
    end
    function randomize!(path::TtPath)
        randomize!(path.agepoints, path.Tpoints, path.boundary, path.detail)
        randomize!(path.constraint, path.boundary)
        randomize!(path.boundary)
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
    
    function initialproposal!(path::TtPath, npoints::Int, dTmax::Number; nattempts = 1_000_000) 
        for _ in 1:nattempts
            randomize!(path)
            collectaccepted!(path, npoints)
            if maxdiff(path.Tsteps) < dTmax
                break
            end
        end
        if maxdiff(path.Tsteps) > dTmax
            @warn "Could not generate initial proposal to satisfy max reheating rate in $nattempts attempts"
        end
        return path
    end

    # Adjust kinetic models
    function movekinetics(zdm::ZRDAAM)
        p = 0.5
        ZRDAAM(
            DzEa = (rand()<p) ? exp(log(zdm.DzEa)+randn()*zdm.DzEa_logsigma) : zdm.DzEa,
            DzD0 = (rand()<p) ? exp(log(zdm.DzD0)+randn()*zdm.DzD0_logsigma) : zdm.DzD0,
            DN17Ea = (rand()<p) ? exp(log(zdm.DN17Ea)+randn()*zdm.DN17Ea_logsigma) : zdm.DN17Ea,
            DN17D0 = (rand()<p) ? exp(log(zdm.DN17D0)+randn()*zdm.DN17D0_logsigma) : zdm.DN17D0,
            rmr0 = (rand()<p) ? reflecting(zdm.rmr0 + randn()*zdm.rmr0_sigma, 0, 1) : zdm.rmr0,
        )
    end
    function movekinetics(adm::RDAAM)
        p = 0.5
        rmr0 = (rand()<p) ? reflecting(adm.rmr0 + randn()*adm.rmr0_sigma, 0, 1) : adm.rmr0
        RDAAM(
            D0L = (rand()<p) ? exp(log(adm.D0L)+randn()*adm.D0L_logsigma) : adm.D0L,
            EaL = (rand()<p) ? exp(log(adm.EaL)+randn()*adm.EaL_logsigma) : adm.EaL,
            EaTrap = (rand()<p) ? exp(log(adm.EaTrap)+randn()*adm.EaTrap_logsigma) : adm.EaTrap,
            rmr0 = rmr0,
            kappa = adm.kappa_rmr0 - rmr0,
        )
    end

    # Adjust model uncertainties of chronometers
    function movesigma!(σcalc::AbstractVector{T}, chrons::AbstractVector{<:Chronometer}) where {T<:AbstractFloat}
        for C in (ZirconFT, ApatiteFT, ZirconHe, ApatiteHe, GenericHe, GenericAr,)
            r = abs(randn(T))
            for i in eachindex(σcalc, chrons)
                if chrons[i] isa C
                    σcalc[i] *= r
                end
            end
        end
        return σcalc
    end

    # Utility function to calculate model ages for all chronometers at once
    function model!(μcalc::AbstractVector{T}, σcalc::AbstractVector{T}, data::Vector{<:Chronometer{T}}, Tsteps::AbstractVector{T}, zdm::ZirconHeliumModel{T}, adm::ApatiteHeliumModel{T}, zftm::ZirconAnnealingModel{T}, mftm::MonaziteAnnealingModel{T}, aftm::ApatiteAnnealingModel{T}; trackhist::Bool=false) where {T<:AbstractFloat}
        imax = argmax(i->length(data[i].agesteps), eachindex(data))
        tsteps = data[imax].tsteps
        tmax = last(tsteps)
        dt = step(tsteps)
        @assert issorted(tsteps)
        @assert eachindex(tsteps) == eachindex(Tsteps)

        # Pre-anneal ZRDAAM samples, if any
        isa(zdm, ZRDAAM) && anneal!(data, ZirconHe{T}, tsteps, Tsteps, zdm)
        # Pre-anneal RDAAM samples, if any
        isa(adm, RDAAM) && anneal!(data, ApatiteHe{T}, tsteps, Tsteps, adm)
        
        # Cycle through each Chronometer, model and calculate log likelihood
        ll = zero(T)
        for i in eachindex(data, μcalc, σcalc)
            c = data[i]
            first_index = 1 + Int((tmax - last(c.tsteps))÷dt)
            if isa(c, GenericAr) || isa(c, GenericHe)
                μcalc[i] = modelage(c, @views(Tsteps[first_index:end]))
                ll += norm_ll(μcalc[i], σcalc[i], val(c), err(c))
            elseif isa(c, ZirconHe)
                μcalc[i] = modelage(c, @views(Tsteps[first_index:end]), zdm)
                ll += norm_ll(μcalc[i], σcalc[i], val(c), err(c))
            elseif isa(c, ApatiteHe)
                μcalc[i] = modelage(c, @views(Tsteps[first_index:end]), adm)
                ll += norm_ll(μcalc[i], σcalc[i], val(c), err(c))
            elseif isa(c, ZirconFT)
                μcalc[i] = modelage(c, @views(Tsteps[first_index:end]), zftm)
                ll += norm_ll(μcalc[i], σcalc[i], val(c), err(c))
            elseif isa(c, MonaziteFT)
                μcalc[i] = modelage(c, @views(Tsteps[first_index:end]), mftm)
                ll += norm_ll(μcalc[i], σcalc[i], val(c), err(c))
            elseif isa(c, ApatiteFT)
                μcalc[i] = modelage(c, @views(Tsteps[first_index:end]), aftm)
                ll += norm_ll(μcalc[i], σcalc[i], val(c), err(c))
            elseif isa(c, ApatiteTrackLength)
                μcalc[i], _ = modellength(c, @views(Tsteps[first_index:end]), aftm; trackhist)
                ll += model_ll(c)
            else
                # NaN if not calculated
                μcalc[i] = T(NaN)
            end
        end

        if trackhist
            # Additional log likelihood term from comparing observed and expected fission track length histograms
            ll += tracklength_histogram_ll!(data, ApatiteTrackLength)
        end

        return ll
    end

    function tracklength_histogram_ll!(data::Vector{<:Chronometer{T}}, ::Type{C}) where {T<:AbstractFloat, C<:ApatiteTrackLength}
        # Initial log likelihood
        ll = zero(T)
        # See how many tracks of type C we have
        n_tracks = count(x->isa(x, C), data)
        if n_tracks > 1
            # Predicted track counts
            i1 = findfirst(x->isa(x,C), data)
            predicted = (data[i1]::C).ldist::Vector{T}
            for i in eachindex(data)
                if isa(data[i], C) && (i > i1)
                    c = data[i]::C
                    predicted .+= c.ldist
                end
            end
            # Scale to match histogram of observed track counts, if any
            Σpredicted = sum(predicted)
            (Σpredicted > 0) && (predicted .*= n_tracks/Σpredicted)

            # Observed track counts
            i2 = findnext(x->isa(x,C), data, i1+1)
            binedges = (data[i2]::C).ledges::FloatRange
            observed = (data[i2]::C).ldist::Vector{T}
            fill!(observed, zero(T))
            for i in eachindex(data)
                if isa(data[i], C) 
                    c = data[i]::C
                    li = Int((val(c) - first(binedges)) ÷ step(binedges)) + firstindex(observed)
                    if firstindex(observed) <= li <= lastindex(observed)
                        observed[li] += 1
                    end
                end
            end

            # Compare observed and predicted length histograms via counting statistics of each bin
            # Approximately equivalent to a Pearson's chi-squared test on the two distributions
            for i in eachindex(observed, predicted)
                # Arrival distributions in each length bin should follow Possson distributions
                ll += logcdf(Poisson(predicted[i]), observed[i])
            end
        end
        return ll
    end
    
## --- Ensure non-allocation of linear algebra
function lu!(A::Tridiagonal{T,V}, pivot::Union{RowMaximum,NoPivot} = RowMaximum();
        check::Bool = true, allowsingular::Bool = false) where {T,V}
    n = size(A, 1)
    has_du2_defined = isdefined(A, :du2) && length(A.du2) == max(0, n-2)
    if has_du2_defined
        du2 = A.du2::V
    else
        du2 = similar(A.d, max(0, n-2))::V
    end
    _lu_tridiag!(Tridiagonal{T,V}(A.dl, A.d, A.du, du2), Vector{BlasInt}(undef, n), pivot, check, allowsingular)
end
function lu!(F::LU{<:Any,<:Tridiagonal}, A::Tridiagonal, pivot::Union{RowMaximum,NoPivot} = RowMaximum();
        check::Bool = true, allowsingular::Bool = false)
    B = F.factors
    size(B) == size(A) || throw(DimensionMismatch())
    copyto!(B, A)
    _lu_tridiag!(B, F.ipiv, pivot, check, allowsingular)
end
@inline function _lu_tridiag!(B::Tridiagonal{T,V}, ipiv, pivot, check, allowsingular) where {T,V}

    # Extract values
    dl = B.dl::V
    d = B.d::V
    du = B.du::V
    du2 = B.du2::V
    n = length(d)

    # Initialize variables
    info = 0
    fill!(du2, 0)

    @inbounds begin
        for i = 1:n
            ipiv[i] = i
        end
        for i = 1:n-2
            # pivot or not?
            if pivot === NoPivot() || abs(d[i]) >= abs(dl[i])
                # No interchange
                if d[i] != 0
                    fact = dl[i]/d[i]
                    dl[i] = fact
                    d[i+1] -= fact*du[i]
                    du2[i] = 0
                end
            else
                # Interchange
                fact = d[i]/dl[i]
                d[i] = dl[i]
                dl[i] = fact
                tmp = du[i]
                du[i] = d[i+1]
                d[i+1] = tmp - fact*d[i+1]
                du2[i] = du[i+1]
                du[i+1] = -fact*du[i+1]
                ipiv[i] = i+1
            end
        end
        if n > 1
            i = n-1
            if pivot === NoPivot() || abs(d[i]) >= abs(dl[i])
                if d[i] != 0
                    fact = dl[i]/d[i]
                    dl[i] = fact
                    d[i+1] -= fact*du[i]
                end
            else
                fact = d[i]/dl[i]
                d[i] = dl[i]
                dl[i] = fact
                tmp = du[i]
                du[i] = d[i+1]
                d[i+1] = tmp - fact*d[i+1]
                ipiv[i] = i+1
            end
        end
        # check for a zero on the diagonal of U
        for i = 1:n
            if d[i] == 0
                info = i
                break
            end
        end
    end
    check && LinearAlgebra._check_lu_success(info, allowsingular)
    return LU{T,Tridiagonal{T,V},typeof(ipiv)}(B, ipiv, convert(LinearAlgebra.BlasInt, info))
end
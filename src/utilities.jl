## ---  Various useful internal utility functions

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

            # Let omega be is the solid angle of intersection normalized by 4pi,
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
        fintlast = intersectionfraction(first(redges),ralpha,d)
        @inbounds for i ∈ I
            # Integrated intersection fraction for each concentric sphere (redges) of crystal
            fint = intersectionfraction(redges[i+1],ralpha,d)

            # Intersection fraction for each spherical shell of the crystal (subtracting
            # one concentric sphere from the next) normalized by shell volume
            dint[i] = (fint-fintlast) / rvolumes[i]
            fintlast = fint
        end
        return dint
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

    function logmeantemp(Tsteps::AbstractArray{T}) where T
        Tf = float(T)
        lΣ = zero(Tf)
        @inbounds for i in eachindex(Tsteps)
            lΣ += log(Tsteps[i]+ Tf(273.15))
        end
        return exp(lΣ/length(Tsteps)) - Tf(273.15)
    end

    """
    ```julia
    simannealsigma(n::Integer, σₑ::Number, σₐ::Number, λₐ::Number)
    ```
    To avoid getting stuck in local optima, combines empirically observed 
    uncertainty `σₑ` in quadrature with an annealling uncertainty which 
    slowly declines from `σₐ` to `0` with a decay constant of `λₐ`, or in  
    other words:

        σannealed = sqrt(σₑ^2 + (σₐ*exp(-λ*n))^2)

    """
    function simannealsigma(n::Integer, σₑ::Number, σₐ::Number, λₐ::Number)
        @assert λₐ >= 0
        sqrt(σₑ^2 + (σₐ * exp(-λₐ*n))^2)
    end
    export simannealsigma


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

    function boundtime(t::Number, boundary::Boundary)
        tmin, tmax = extrema(boundary.agepoints)
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
        Tmin, Tmax = extrema(boundary.T₀)
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
    function movepoint!(agepointsₚ::Vector{T}, Tpointsₚ::Vector{T}, k::Int, σⱼt::T, σⱼT::T, boundary::Boundary{T}) where {T}

        # Move the age of one model point
        agepointsₚ[k] += randn() * σⱼt

        # Move the Temperature of one model point
        Tpointsₚ[k] += randn() * σⱼT

        # Apply time boundary conditions
        agepointsₚ[k] = boundtime(agepointsₚ[k], boundary)

        # Apply Temperature boundary conditions
        Tpointsₚ[k] = boundtemp(Tpointsₚ[k], boundary)

        return agepointsₚ[k], Tpointsₚ[k]
    end

    function addpoint!(agepointsₚ::Vector{T}, Tpointsₚ::Vector{T}, σⱼt::Vector{T}, σⱼT::Vector{T}, k::Int, boundary::Boundary{T}) where {T}
        @assert eachindex(agepointsₚ) == eachindex(Tpointsₚ) == eachindex(σⱼt) == eachindex(σⱼT)

        tmin, tmax = extrema(boundary.agepoints)
        Tmin, Tmax = extrema(boundary.T₀)

        # Pick an age uniformly within the boundaries
        agepointsₚ[k] = tmin + rand()*(tmax-tmin)

        # Find the closest existing points (if any)
        ages = view(agepointsₚ, 1:k-1)
        i₋ = findclosestbelow(agepointsₚ[k], ages)
        i₊ = findclosestabove(agepointsₚ[k], ages)

        # Find values for the closest younger point
        inbounds₋ = firstindex(ages) <= i₋ <= lastindex(ages)
        t₋ = inbounds₋ ? agepointsₚ[i₋] : tmin
        T₋ = inbounds₋ ? Tpointsₚ[i₋] : Tmin
        σⱼt₋ = inbounds₋ ? σⱼt[i₋] : (tmax-tmin)/60
        σⱼT₋ = inbounds₋ ? σⱼT[i₋] : (Tmax-Tmin)/60

        # Find values for the closest older point
        inbounds₊ = firstindex(ages) <= i₊ <= lastindex(ages)
        t₊ = inbounds₊ ? agepointsₚ[i₊] : tmax
        T₊ = inbounds₊ ? Tpointsₚ[i₊] : Tmax
        σⱼt₊ = inbounds₊ ? σⱼt[i₊] : (tmax-tmin)/60
        σⱼT₊ = inbounds₊ ? σⱼT[i₊] : (Tmax-Tmin)/60

        # Interpolate
        f = (agepointsₚ[k] - t₋) / (t₊ - t₋)
        f *= !isnan(f)
        Tpointsₚ[k] = f*T₊ + (1-f)*T₋
        σⱼt[k] = f*σⱼt₊ + (1-f)*σⱼt₋
        σⱼT[k] = f*σⱼT₊ + (1-f)*σⱼT₋

        # Move the point from the interpolated value
        movepoint!(agepointsₚ, Tpointsₚ, k, σⱼt[k], σⱼT[k], boundary)

        return agepointsₚ[k], Tpointsₚ[k]
    end

    function movebounds!(boundary::Boundary)
        @inbounds for i in eachindex(boundary.Tpointsₚ)
            boundary.Tpointsₚ[i] = boundary.T₀[i] + rand()*boundary.ΔT[i]
        end
        boundary
    end
    function movebounds!(constraint::Constraint, boundary::Boundary)
        @inbounds for i in eachindex(constraint.Tpointsₚ)
            constraint.agepointsₚ[i] = boundtime(rand(constraint.agedist[i]), boundary)
            constraint.Tpointsₚ[i] = boundtemp(rand(constraint.Tdist[i]), boundary)
        end
        constraint
    end

    function movekinetics(zdm::ZRDAAM)
        rn = rand(1:5)
        zdmₚ = ZRDAAM(
            DzEa = (rn==1) ? exp(log(zdm.DzEa)+randn()*zdm.DzEa_logsigma) : zdm.DzEa,
            DzD0 = (rn==2) ? exp(log(zdm.DzD0)+randn()*zdm.DzD0_logsigma) : zdm.DzD0,
            DN17Ea = (rn==3) ? exp(log(zdm.DN17Ea)+randn()*zdm.DN17Ea_logsigma) : zdm.DN17Ea,
            DN17D0 = (rn==4) ? exp(log(zdm.DN17D0)+randn()*zdm.DN17D0_logsigma) : zdm.DN17D0,
            rmr0 = (rn==5) ? reflecting(zdm.rmr0 + randn()*zdm.rmr0_sigma, 0, 1) : zdm.rmr0,
        )
    end
    function movekinetics(adm::RDAAM)
        rn = rand(1:4)
        rmr0 = (rn==4) ? reflecting(adm.rmr0 + randn()*adm.rmr0_sigma, 0, 1) : adm.rmr0
        RDAAM(
            D0L = (rn==1) ? exp(log(adm.D0L)+randn()*adm.D0L_logsigma) : adm.D0L,
            EaL = (rn==2) ? exp(log(adm.EaL)+randn()*adm.EaL_logsigma) : adm.EaL,
            EaTrap = (rn==3) ? exp(log(adm.EaTrap)+randn()*adm.EaTrap_logsigma) : adm.EaTrap,
            rmr0 = rmr0,
            kappa = adm.kappa_rmr0 - rmr0,
        )
    end

    # Normal log likelihood
    @inline function norm_ll(mu::Number, sigma::Number, x::Number)
        δ = x - mu
        σ² = sigma^2
        -0.5*(log(2*pi*σ²) + δ^2/σ²)
    end
    function norm_ll(mu::AbstractVector{T}, sigma::AbstractVector{T}, x::AbstractVector{T}) where {T<:Number}
        ll = zero(float(T))
        for i in eachindex(mu, sigma, x)
            ll += norm_ll(mu[i], sigma[i], x[i])
        end
        return ll
    end
    @inline function norm_ll(mu::Number, sigma::Number, x::Number, x_sigma::Number)
        δ = x - mu
        σ² = sigma^2 + x_sigma^2
        -0.5*(log(2*pi*σ²) + δ^2/σ²)
    end
    function norm_ll(mu::AbstractVector{T}, sigma::AbstractVector{T}, x::AbstractVector{T}, x_sigma::AbstractVector{T}) where {T<:Number}
        ll = zero(float(T))
        for i in eachindex(mu, sigma, x, x_sigma)
            ll += norm_ll(mu[i], sigma[i], x[i], x_sigma[i])
        end
        return ll
    end

    # Specialized ll functions for ZRDAAM and RDAAM
    function loglikelihood(zdmₚ::ZRDAAM, zdm::ZRDAAM)
        norm_ll(log(zdm.DzD0), zdm.DzD0_logsigma, log(zdmₚ.DzD0)) + 
        norm_ll(log(zdm.DzEa), zdm.DzEa_logsigma, log(zdmₚ.DzEa)) + 
        norm_ll(log(zdm.DN17D0), zdm.DN17D0_logsigma, log(zdmₚ.DN17D0)) + 
        norm_ll(log(zdm.DN17Ea), zdm.DN17Ea_logsigma, log(zdmₚ.DN17Ea)) +
        norm_ll(zdm.rmr0, zdm.rmr0_sigma, zdmₚ.rmr0)
    end
    function loglikelihood(admₚ::RDAAM, adm::RDAAM)
        norm_ll(log(adm.D0L), adm.D0L_logsigma, log(admₚ.D0L)) + 
        norm_ll(log(adm.EaL), adm.EaL_logsigma, log(admₚ.EaL)) + 
        norm_ll(log(adm.EaTrap), adm.EaTrap_logsigma, log(admₚ.EaTrap))+
        norm_ll(adm.rmr0, adm.rmr0_sigma, admₚ.rmr0)
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
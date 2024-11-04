## ---  Various useful internal utility functions

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
    function maxdiff(x::AbstractVector{T}) where {T}
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
    function mindiff(x::AbstractVector{T}) where {T}
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
    function maxabsdiff(x::AbstractVector{T}) where {T}
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

    function diff_ll(x::AbstractVector, μ::Number, σ::Number)
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
    function isdistinct(points::DenseArray, npoints::Int, k::Int, δ::Number)
        @inbounds for i = 1:npoints
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
                n += isdistinct(points, npoints, i, δ)
            end
        end
        return n
    end


    """
    ```julia
    simannealsigma(n, σAnalytical; [simannealmodel::NamedTuple])
    simannealsigma(n::Integer, σAnalytical::Number, σmodel::Number, σannealing::Number, λannealing::Number)
    ```
    To avoid getting stuck in local optima, decrease uncertainty slowly by
    simulated annealing. Parameters are specified as a tuple `simannealmodel` of the
    form (σₘ, σᵢ, λ), where annealing uncertainty declines from `σᵢ+σₘ` to `σₘ`
    with a decay constant of λ.

    Returns the annealing uncertainty added in quadrature with analytical
    uncertainty, or in other words

        sigma = sqrt(σAnalytical^2 + (σᵢ*exp(-λ*n) + σₘ)^2)

    """
    function simannealsigma(n::Integer, σAnalytical::Number; simannealmodel::NamedTuple=(σmodel=25.0, σannealing=35.0, λannealing=10/10^5))
        simannealsigma(n, σAnalytical, simannealmodel.σmodel, simannealmodel.σannealing, simannealmodel.λannealing)
    end
    function simannealsigma(n::Integer, σAnalytical::Number, σmodel::Number, σannealing::Number, λannealing::Number)
        σCombined = σannealing * exp(-λannealing*n) + σmodel
        return sqrt(σAnalytical^2 + σCombined^2)
    end
    export simannealsigma


    # Move a t-T point and apply boundary conditions
    function movepoint!(agepointsₚ, Tpointsₚ, k, tmin, tmax, Tmin, Tmax, σⱼt, σⱼT, boundarytype::Symbol=:hard)
        # Ensure bounds are in proper order
        tmin < tmax || ((tmin, tmax) = (tmax, tmin))
        Tmin < Tmax || ((Tmin, Tmax) = (Tmax, Tmin))

        # Move the age of one model point
        agepointsₚ[k] += randn() * σⱼt

        # Move the Temperature of one model point
        Tpointsₚ[k] += randn() * σⱼT

        # Reflecting boundary conditions
        if boundarytype === :reflecting
            if agepointsₚ[k] < tmin
                agepointsₚ[k] = tmin - (agepointsₚ[k] - tmin)
            elseif agepointsₚ[k] > tmax
                agepointsₚ[k] = tmax - (agepointsₚ[k] - tmax)
            end

            if Tpointsₚ[k] < Tmin
                Tpointsₚ[k] = Tmin - (Tpointsₚ[k] - Tmin)
            elseif Tpointsₚ[k] > Tmax
                Tpointsₚ[k] = Tmax - (Tpointsₚ[k] - Tmax)
            end
        end

        # Hard boundary conditions
        agepointsₚ[k] = min(max(agepointsₚ[k], tmin), tmax)
        Tpointsₚ[k] = min(max(Tpointsₚ[k], Tmin),Tmax)

        return  agepointsₚ, Tpointsₚ
    end

    function movekinetics(zdm::ZRDAAM)
        rn = rand(1:4)
        zdmₚ = ZRDAAM(
            DzEa = (rn==1) ? exp(log(zdm.DzEa)+randn()*zdm.DzEa_logsigma) : zdm.DzEa,
            DzD0 = (rn==2) ? exp(log(zdm.DzD0)+randn()*zdm.DzD0_logsigma) : zdm.DzD0,
            DN17Ea = (rn==3) ? exp(log(zdm.DN17Ea)+randn()*zdm.DN17Ea_logsigma) : zdm.DN17Ea,
            DN17D0 = (rn==4) ? exp(log(zdm.DN17D0)+randn()*zdm.DN17D0_logsigma) : zdm.DN17D0,
        )
    end
    function movekinetics(adm::RDAAM)
        rn = rand(1:3)
        RDAAM(
            D0L = (rn==1) ? exp(log(adm.D0L)+randn()*adm.D0L_logsigma) : adm.D0L,
            EaL = (rn==2) ? exp(log(adm.EaL)+randn()*adm.EaL_logsigma) : adm.EaL,
            EaTrap = (rn==3) ? exp(log(adm.EaTrap)+randn()*adm.EaTrap_logsigma) : adm.EaTrap,
        )
    end
    
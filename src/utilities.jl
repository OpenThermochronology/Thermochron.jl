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

    # Move a t-T point and apply boundary conditions
    function movepoint!(agepointsₚ::Vector{T}, Tpointsₚ::Vector{T}, k::Int, σⱼt::T, σⱼT::T, boundary::Boundary{T}) where {T}

        # Move the age of one model point
        agepointsₚ[k] += randn() * σⱼt

        # Move the Temperature of one model point
        Tpointsₚ[k] += randn() * σⱼT

        # Apply time boundary conditions
        tmin, tmax = extrema(boundary.agepoints)
        if boundary.tboundary === :reflecting
            agepointsₚ[k] = reflecting(agepointsₚ[k], tmin, tmax)
        elseif boundary.tboundary === :hard
            agepointsₚ[k] = hard(agepointsₚ[k], tmin, tmax)
        elseif boundary.tboundary === :periodic
            agepointsₚ[k] = periodic(agepointsₚ[k], tmin, tmax)
        else
            @error "`tboundary` must be either `:reflecting`, `:hard`, or `:periodic`"
        end

        # Apply Temperature boundary conditions
        Tmin, Tmax = extrema(boundary.T₀)
        if boundary.Tboundary === :reflecting
            Tpointsₚ[k] = reflecting(Tpointsₚ[k], Tmin, Tmax)
        elseif boundary.Tboundary === :hard
            Tpointsₚ[k] = hard(Tpointsₚ[k], Tmin, Tmax)
        elseif boundary.Tboundary === :periodic
            Tpointsₚ[k] = periodic(Tpointsₚ[k], Tmin, Tmax)
        else
            @error "`Tboundary` must be either `:reflecting`, `:hard`, or `:periodic`"
        end

        return agepointsₚ[k], Tpointsₚ[k]
    end

    function movebounds!(boundary::Boundary)
        @inbounds for i in eachindex(boundary.Tpointsₚ)
            boundary.Tpointsₚ[i] = boundary.T₀[i] + rand()*boundary.ΔT[i]
        end
        boundary
    end
    function movebounds!(constraint::Constraint)
        @inbounds for i in eachindex(constraint.Tpointsₚ)
            constraint.agepointsₚ[i] = constraint.Age₀[i] + rand()*constraint.ΔAge[i]
            constraint.Tpointsₚ[i] = constraint.T₀[i] + rand()*constraint.ΔT[i]
        end
        constraint
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
    
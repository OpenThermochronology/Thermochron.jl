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

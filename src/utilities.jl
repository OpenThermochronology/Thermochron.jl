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

## --- Modify LinearAlgebra.lu! to reuse du2 for pivoting

# using LinearAlgebra: BlasInt, checknonsingular
# # Modified from LinearAlgebra stdlib to reuse `du2` and `ipiv`
# function lu!(A::Tridiagonal{T,V}, ipiv::Vector{BlasInt}, pivot::Union{RowMaximum,NoPivot} = RowMaximum(); check::Bool = true) where {T,V}

#     # Extract values
#     n = size(A, 1)
#     @assert length(ipiv) == n

#     # initialize variables
#     info = 0
#     dl = A.dl
#     d = A.d
#     du = A.du
#     if dl === du
#         throw(ArgumentError("off-diagonals of `A` must not alias"))
#     end
#     # Check if Tridiagonal matrix already has du2 for pivoting
#     d2len = max(0, n-2) # Proper length of a second off-diagonal
#     has_du2_defined = isdefined(A, :du2) && length(A.du2) == d2len
#     du2 = (has_du2_defined ? A.du2 : similar(d, d2len))::V
#     fill!(du2, 0)

#     @inbounds begin
#         for i = 1:n
#             ipiv[i] = i
#         end
#         for i = 1:n-2
#             # pivot or not?
#             if pivot === NoPivot() || abs(d[i]) >= abs(dl[i])
#                 # No interchange
#                 if d[i] != 0
#                     fact = dl[i]/d[i]
#                     dl[i] = fact
#                     d[i+1] -= fact*du[i]
#                     du2[i] = 0
#                 end
#             else
#                 # Interchange
#                 fact = d[i]/dl[i]
#                 d[i] = dl[i]
#                 dl[i] = fact
#                 tmp = du[i]
#                 du[i] = d[i+1]
#                 d[i+1] = tmp - fact*d[i+1]
#                 du2[i] = du[i+1]
#                 du[i+1] = -fact*du[i+1]
#                 ipiv[i] = i+1
#             end
#         end
#         if n > 1
#             i = n-1
#             if pivot === NoPivot() || abs(d[i]) >= abs(dl[i])
#                 if d[i] != 0
#                     fact = dl[i]/d[i]
#                     dl[i] = fact
#                     d[i+1] -= fact*du[i]
#                 end
#             else
#                 fact = d[i]/dl[i]
#                 d[i] = dl[i]
#                 dl[i] = fact
#                 tmp = du[i]
#                 du[i] = d[i+1]
#                 d[i+1] = tmp - fact*d[i+1]
#                 ipiv[i] = i+1
#             end
#         end
#         # check for a zero on the diagonal of U
#         for i = 1:n
#             if d[i] == 0
#                 info = i
#                 break
#             end
#         end
#     end
#     B = has_du2_defined ? A : Tridiagonal{T,V}(dl, d, du, du2)
#     check && checknonsingular(info, pivot)
#     return LU{T,Tridiagonal{T,V},typeof(ipiv)}(B, ipiv, convert(BlasInt, info))
# end

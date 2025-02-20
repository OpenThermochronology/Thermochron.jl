
struct MultipleDomain{T<:AbstractFloat, C<:Union{SphericalAr{T}, PlanarAr{T}}} <: AbsoluteChronometer{T}
    age::Vector{T}
    age_sigma::Vector{T}
    fraction_released::Vector{T}
    tsteps_experimental::Vector{T}
    Tsteps_experimental::Vector{T}
    fit::BitVector
    offset::T
    domains::Vector{C}
    volume_fraction::Vector{T}
    model_age::Vector{T}
    model_parent::Vector{T}
    model_daughter::Vector{T}
    tsteps_degassing::FloatRange
    Tsteps_degassing::Vector{T}
    agesteps::FloatRange
end
function MultipleDomain(T=Float64, C=PlanarAr;
        age::AbstractVector,
        age_sigma::AbstractVector,
        fraction_released::AbstractVector,
        tsteps_experimental::AbstractVector,
        Tsteps_experimental::AbstractVector,
        fit::AbstractVector,
        offset::Number = zero(T),
        Ea::AbstractVector,
        lnD0_a_2::AbstractVector,
        volume_fraction::AbstractVector,
        r::Number = 100,
        dr::Number = one(T),
        K40::Number = 16.34, 
        agesteps::AbstractRange,
    )
    # Check input arrays are the right size and ordered properly
    @assert eachindex(Ea) == eachindex(lnD0_a_2) == eachindex(volume_fraction)
    @assert eachindex(age) == eachindex(age_sigma) == eachindex(fraction_released) == eachindex(tsteps_experimental) == eachindex(Tsteps_experimental)
    @assert issorted(tsteps_experimental, lt=<=) "Degassing time steps must be in strictly increasing order"

    # Interpolate degassing t-T steps to the same resolution as the 
    tsteps_degassing = range(Float64(first(tsteps_experimental)), Float64(last(tsteps_experimental)), length=length(agesteps))
    Tsteps_degassing = linterp1(tsteps_experimental, T.(Tsteps_experimental), tsteps_degassing) 
    model_age = zeros(T, length(tsteps_degassing))
    model_parent = zeros(T, length(tsteps_degassing))
    model_daughter = zeros(T, length(tsteps_degassing))

    # Ensure volume fraction sums to one
    volume_fraction ./= nansum(volume_fraction) 
    
    # Allocate domains
    D0 = @. exp(lnD0_a_2)*(r/10000)^2
    domains = [C(T; age=nanmean(age), age_sigma=nanstd(age), offset, r, dr, D0=D0[i], Ea=Ea[i], K40, agesteps) for i in eachindex(D0, Ea)]
    return MultipleDomain{T,C{T}}(
        T.(age),
        T.(age_sigma),
        T.(fraction_released),
        T.(tsteps_experimental),
        T.(Tsteps_experimental),
        Bool.(fit),
        T(offset),
        domains,
        T.(volume_fraction),
        model_age,
        model_parent,
        model_daughter,
        tsteps_degassing,
        Tsteps_degassing,
        floatrange(agesteps),
    )
end

function modelage(mdd::MultipleDomain{T}, Tsteps::AbstractVector; redegasparent::Bool=false) where {T<:AbstractFloat}
    age = fill!(mdd.model_age, zero(T))
    parent = fill!(mdd.model_parent, zero(T))
    daughter = fill!(mdd.model_daughter, zero(T))
    for i in eachindex(mdd.domains, mdd.volume_fraction)
        domain = mdd.domains[i]
        modelage(domain, Tsteps)
        p, d = degas!(domain, mdd.tsteps_degassing, mdd.Tsteps_degassing; redegasparent)
        @. parent += p * mdd.volume_fraction[i]
        @. daughter += d * mdd.volume_fraction[i]
    end
    for i in eachindex(parent, daughter)
        age[i] = newton_ar_age(daughter[i], parent[i])
    end
    parent ./= nansum(parent)
    nancumsum!(parent)
    return age, parent
end

function degas!(mineral::PlanarAr{T}, tsteps_degassing::FloatRange, Tsteps_degassing::AbstractVector{T}; redegasparent::Bool=false) where T <: AbstractFloat

    # Damage and annealing constants
    D0 = (mineral.D0*10000^2)::T              # cm^2/sec, converted to micron^2/sec  
    Ea = mineral.Ea::T                      # kJ/mol
    R = 0.008314472                         # kJ/(K*mol)
    ΔT = mineral.offset::T + 273.15         # Conversion from C to K, plus temperature offset from the

    # Diffusivities of crystalline and amorphous endmembers
    De = mineral.De::Vector{T}
    @assert firstindex(De) == firstindex(Tsteps_degassing)
    @assert lastindex(De) >= lastindex(Tsteps_degassing)
    @turbo for i ∈ eachindex(Tsteps_degassing)
        De[i] = D0 * exp(-Ea / R / (Tsteps_degassing[i] + ΔT)) # micron^2/Myr
    end

    # Get time and radius discretization
    dr = step(mineral.rsteps)
    nrsteps = mineral.nrsteps::Int
    dt = step(tsteps_degassing)
    ntsteps = length(tsteps_degassing)
    step_parent = @views(mineral.step_parent[1:ntsteps])
    step_daughter = @views(mineral.step_daughter[1:ntsteps])
    @assert eachindex(tsteps_degassing) == eachindex(Tsteps_degassing) == Base.OneTo(ntsteps)

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}
    fill!(β, 2 * dr^2 / (De[1]*dt))

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = mineral.y

    # Output matrix for all timesteps
    # No coordinate transform required for slab geometry, so here u is the diffusing Ar profile
    u = mineral.u::DenseMatrix{T}
    fill!(u, zero(T)) 
    u[:,1] = y # Start with final Ar profile from end of last geologic inversion

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = mineral.A       # Tridiagonal matrix
    F = mineral.F       # LU object for in-place lu factorization
    
    daughterᵢ₋ = nanmean(@views(u[2:end-1, 1]))
    @inbounds for i in Base.OneTo(ntsteps)

        # Update β for current temperature
        fill!(β, 2 * dr^2 / (De[i]*dt))

        # Update tridiagonal matrix
        fill!(A.dl, 1)         # Sub-diagonal
        @. A.d = -2 - β        # Diagonal
        fill!(A.du, 1)         # Supra-diagonal

        # Neumann inner boundary condition (-u(i,1) + u(i,2) = 0)
        A.du[1] = -1
        A.d[1] = 1
        y[1] = 0

        # Dirichlet outer boundary condition (u(i,end) = u(i-1,end))
        A.dl[nrsteps-1] = 0
        A.d[nrsteps] = 1
        y[nrsteps] = u[nrsteps,i]

        # RHS of tridiagonal Crank-Nicholson equation for regular grid cells.
        # From Ketcham, 2005 https://doi.org/10.2138/rmg.2005.58.11
        @turbo for k = 2:nrsteps-1
            𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i+1] = y

        daughterᵢ = nanmean(@views(u[2:end-1, i+1]))
        step_daughter[i] = daughterᵢ₋ - daughterᵢ
        daughterᵢ₋ = daughterᵢ
    end

    # Now diffuse parent isotope tracer, (as Ar-39), if neccesary
    if redegasparent || !(0 < sum(step_parent))

        # Convert diffusivity from that of Ar-40 to that of Ar-39 given Dₗ/Dₕ ~ (mₕ/mₗ)^β 
        # C.f. Luo et al. 2021 (doi: 10.7185/geochemlet.2128) for He in albite melt
        # β=0.355 ± 0.012 at 3000 K, decreasing to β=0.322 ± 0.019 at 1700 K
        De .*= (40/39)^0.3
        
        # Initialize u matrix
        fill!(u, zero(T)) 
        u[2:end-1,1] = mineral.r40K
        u[1,1] = u[2,1]
        u[end,1] = 0

        parentᵢ₋ = nanmean(@views(u[2:end-1, 1]))
        @inbounds for i = 1:ntsteps

            # Update β for current temperature
            fill!(β, 2 * dr^2 / (De[i]*dt))

            # Update tridiagonal matrix
            fill!(A.dl, 1)         # Sub-diagonal
            @. A.d = -2 - β        # Diagonal
            fill!(A.du, 1)         # Supra-diagonal

            # Neumann inner boundary condition (-u(i,1) + u(i,2) = 0)
            A.du[1] = -1
            A.d[1] = 1
            y[1] = 0

            # Dirichlet outer boundary condition (u(i,end) = u(i-1,end))
            A.dl[nrsteps-1] = 0
            A.d[nrsteps] = 1
            y[nrsteps] = u[nrsteps,i]

            # RHS of tridiagonal Crank-Nicholson equation for regular grid cells.
            # From Ketcham, 2005 https://doi.org/10.2138/rmg.2005.58.11
            @turbo for k = 2:nrsteps-1
                𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
                y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊
            end

            # Invert using tridiagonal matrix algorithm
            # equivalent to u[:,i] = A\y
            lu!(F, A, allowsingular=true)
            ldiv!(F, y)
            u[:,i+1] = y

            parentᵢ = nanmean(@views(u[2:end-1, i+1]))
            step_parent[i] = parentᵢ₋ - parentᵢ
            parentᵢ₋ = parentᵢ
        end
    end

    # Return views of the resulting step ages and degassing fractions
    return step_parent, step_daughter
end
function degas!(mineral::SphericalAr{T}, tsteps_degassing::FloatRange, Tsteps_degassing::AbstractVector{T}; redegasparent::Bool=false) where T <: AbstractFloat

    # Damage and annealing constants
    D0 = (mineral.D0*10000^2)::T              # cm^2/sec, converted to micron^2/sec  
    Ea = mineral.Ea::T                      # kJ/mol
    R = 0.008314472                         # kJ/(K*mol)
    ΔT = mineral.offset::T + 273.15         # Conversion from C to K, plus temperature offset from the

    # Diffusivities of crystalline and amorphous endmembers
    De = mineral.De::Vector{T}
    @assert firstindex(De) == firstindex(Tsteps_degassing)
    @assert lastindex(De) >= lastindex(Tsteps_degassing)
    @turbo for i ∈ eachindex(Tsteps_degassing)
        De[i] = D0 * exp(-Ea / R / (Tsteps_degassing[i] + ΔT)) # micron^2/Myr
    end

    # Get time and radius discretization
    rsteps = mineral.rsteps
    dr = step(mineral.rsteps)
    nrsteps = mineral.nrsteps::Int
    relvolumes = mineral.relvolumes::Vector{T}
    dt = step(tsteps_degassing)
    ntsteps = length(tsteps_degassing)
    step_parent = @views(mineral.step_parent[1:ntsteps])
    step_daughter = @views(mineral.step_daughter[1:ntsteps])
    @assert eachindex(tsteps_degassing) == eachindex(Tsteps_degassing) == Base.OneTo(ntsteps)

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}
    fill!(β, 2 * dr^2 / (De[1]*dt))

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = mineral.y

    # Output matrix for all timesteps
    # No coordinate transform required for slab geometry, so here u is the diffusing Ar profile
    u = mineral.u::DenseMatrix{T}
    fill!(u, zero(T)) 
    u[:,1] = y # Start with final Ar profile from end of last geologic inversion

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = mineral.A       # Tridiagonal matrix
    F = mineral.F       # LU object for in-place lu factorization
    
    @inbounds for i in Base.OneTo(ntsteps)

        # Update β for current temperature
        fill!(β, 2 * dr^2 / (De[i]*dt))

        # Update tridiagonal matrix
        fill!(A.dl, 1)         # Sub-diagonal
        @. A.d = -2 - β        # Diagonal
        fill!(A.du, 1)         # Supra-diagonal

        # Neumann inner boundary condition (u(i,1) + u(i,2) = 0)
        A.du[1] = 1
        A.d[1] = 1
        y[1] = 0

        # Dirichlet outer boundary condition (u(i,end) = u(i-1,end))
        A.dl[nrsteps-1] = 0
        A.d[nrsteps] = 1
        y[nrsteps] = u[nrsteps,i]

        # RHS of tridiagonal Crank-Nicholson equation for regular grid cells.
        # From Ketcham, 2005 https://doi.org/10.2138/rmg.2005.58.11
        @turbo for k = 2:nrsteps-1
            𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i+1] = y
    end

    # Calculate daughter lost to diffusion at each step
    daughterᵢ₋ = nanmean(@views(u[2:end-1, 1])./=rsteps, relvolumes)
    @inbounds for i in Base.OneTo(ntsteps)
        daughterᵢ = nanmean(@views(u[2:end-1, i+1])./=rsteps, relvolumes)
        step_daughter[i] = daughterᵢ₋ - daughterᵢ
        daughterᵢ₋ = daughterᵢ
    end

    # Now diffuse parent isotope tracer, (as Ar-39), if neccesary
    if redegasparent || !(0 < sum(step_parent))
        
        # Convert diffusivity from that of Ar-40 to that of Ar-39 given Dₗ/Dₕ ~ (mₕ/mₗ)^β 
        # C.f. Luo et al. 2021 (doi: 10.7185/geochemlet.2128) for He in albite melt
        # β=0.355 ± 0.012 at 3000 K, decreasing to β=0.322 ± 0.019 at 1700 K
        De .*= (40/39)^0.3

        # Initialize u matrix
        fill!(u, zero(T)) 
        u[2:end-1,1] = mineral.r40K 
        u[2:end-1,1] .*= mineral.rsteps
        u[1,1] = -u[2,1]
        u[end,1] = 0

        @inbounds for i = 1:ntsteps

            # Update β for current temperature
            fill!(β, 2 * dr^2 / (De[i]*dt))

            # Update tridiagonal matrix
            fill!(A.dl, 1)         # Sub-diagonal
            @. A.d = -2 - β        # Diagonal
            fill!(A.du, 1)         # Supra-diagonal

            # Neumann inner boundary condition (u(i,1) + u(i,2) = 0)
            A.du[1] = 1
            A.d[1] = 1
            y[1] = 0

            # Dirichlet outer boundary condition (u(i,end) = u(i-1,end))
            A.dl[nrsteps-1] = 0
            A.d[nrsteps] = 1
            y[nrsteps] = u[nrsteps,i]

            # RHS of tridiagonal Crank-Nicholson equation for regular grid cells.
            # From Ketcham, 2005 https://doi.org/10.2138/rmg.2005.58.11
            @turbo for k = 2:nrsteps-1
                𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
                y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊
            end

            # Invert using tridiagonal matrix algorithm
            # equivalent to u[:,i] = A\y
            lu!(F, A, allowsingular=true)
            ldiv!(F, y)
            u[:,i+1] = y
        end

        # Calculate parent lost to diffusion at each step
        parentᵢ₋ = nanmean(@views(u[2:end-1, 1])./=rsteps, relvolumes)
        @inbounds for i in Base.OneTo(ntsteps)
            parentᵢ = nanmean(@views(u[2:end-1, i+1])./=rsteps, relvolumes)
            step_parent[i] = parentᵢ₋ - parentᵢ
            parentᵢ₋ = parentᵢ
        end
    end

    # Return views of the resulting step ages and degassing fractions
    return step_parent, step_daughter
end
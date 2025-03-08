## -- Multiple domain diffusivity

    """
    ```julia
    MDDiffusivity(
        D0::NTuple{N,T}             # [cm^2/sec] Maximum diffusivity
        D0_logsigma::NTuple{N,T}    # [unitless] log uncertainty (default = 1/2 = a factor of ℯ two-sigma)
        Ea::T                       # [kJ/mol] Activation energy
        Ea_logsigma::T              # [unitless] log uncertainty (default = 1/2 = a factor of ℯ two-sigma)
    )
    ```
    Multiple diffusivities for multiple domains
    """
    Base.@kwdef struct MDDiffusivity{T<:AbstractFloat, N} <: DiffusivityModel{T}
        D0::NTuple{N,T}             # [cm^2/sec] Maximum diffusivity
        D0_logsigma::NTuple{N,T}    # [unitless] log uncertainty (default = 1/2 = a factor of ℯ two-sigma)
        Ea::T                       # [kJ/mol] Activation energy
        Ea_logsigma::T              # [unitless] log uncertainty (default = 1/2 = a factor of ℯ two-sigma)
    end
    Base.getindex(d::MDDiffusivity{T}, i::Int) where {T} = Diffusivity{T}(d.D0[i], d.D0_logsigma[i], d.Ea, d.Ea_logsigma)

    # Query MDDiffusivities from a KineticResult
    function MDDiffusivity(kr::KineticResult)
        any(x->isa(x, MDDiffusivity), kr) || return nothing
        ia = findall(x->isa(x, MDDiffusivity), kr[:,1])
        return collect(kr[ia,:]')
    end
    
## -- Multiple domain diffusion functions

function degas!(mineral::PlanarAr{T}, tsteps_degassing::FloatRange, Tsteps_degassing::AbstractVector{T}, dm::Diffusivity{T}; fuse::Bool=true, redegasparent::Bool=false) where T <: AbstractFloat

    # Damage and annealing constants
    D0 = (dm.D0*10000^2)::T                 # cm^2/sec, converted to micron^2/sec  
    Ea = dm.Ea::T                           # kJ/mol
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
    @inbounds for i in Base.OneTo(ntsteps-fuse)

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
        step_daughter[i] = max(daughterᵢ₋ - daughterᵢ, zero(T))
        daughterᵢ₋ = daughterᵢ
    end
    # Degas all remaining daughter in last step if fuse==true
    fuse && (step_daughter[ntsteps] = daughterᵢ₋)

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
        @inbounds for i in Base.OneTo(ntsteps-fuse)

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
            step_parent[i] = max(parentᵢ₋ - parentᵢ, zero(T))
            parentᵢ₋ = parentᵢ
        end
        # Degas all remaining parent in last step if fuse==true
        fuse && (step_parent[ntsteps] = parentᵢ₋)
    end

    # Return views of the resulting step ages and degassing fractions
    return step_parent, step_daughter
end
function degas!(mineral::SphericalAr{T}, tsteps_degassing::FloatRange, Tsteps_degassing::AbstractVector{T}, dm::Diffusivity{T}; fuse::Bool=true, redegasparent::Bool=false) where T <: AbstractFloat

    # Damage and annealing constants
    D0 = (dm.D0*10000^2)::T                 # cm^2/sec, converted to micron^2/sec  
    Ea = dm.Ea::T                           # kJ/mol
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
    
    @inbounds for i in Base.OneTo(ntsteps-fuse)

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
    @inbounds for i in Base.OneTo(ntsteps-fuse)
        daughterᵢ = nanmean(@views(u[2:end-1, i+1])./=rsteps, relvolumes)
        step_daughter[i] = max(daughterᵢ₋ - daughterᵢ, zero(T))
        daughterᵢ₋ = daughterᵢ
    end
    # Degas all remaining daughter in last step if fuse==true
    fuse && (step_daughter[ntsteps] = daughterᵢ₋)

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

        @inbounds for i in Base.OneTo(ntsteps-fuse)

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
        @inbounds for i in Base.OneTo(ntsteps-fuse)
            parentᵢ = nanmean(@views(u[2:end-1, i+1])./=rsteps, relvolumes)
            step_parent[i] = max(parentᵢ₋ - parentᵢ, zero(T))
            parentᵢ₋ = parentᵢ
        end
        # Degas all remaining parent in last step if fuse==true
        fuse && (step_parent[ntsteps] = parentᵢ₋)
    end

    # Return views of the resulting step ages and degassing fractions
    return step_parent, step_daughter
end


function modelage(mdd::MultipleDomain{T}, Tsteps::AbstractVector, dm::MDDiffusivity{T}; redegasparent::Bool=false) where {T<:AbstractFloat}
    age = fill!(mdd.model_age, zero(T))
    parent = fill!(mdd.model_parent, zero(T))
    daughter = fill!(mdd.model_daughter, zero(T))
    fraction = fill!(mdd.model_fraction, zero(T))
    fuse = mdd.fuse::Bool
    # Degas
    for i in eachindex(mdd.domains, mdd.volume_fraction)
        domain = mdd.domains[i]
        modelage(domain, Tsteps, dm[i])
        p, d = degas!(domain, mdd.tsteps_degassing, mdd.Tsteps_degassing, dm[i]; fuse, redegasparent)
        @. parent += p * mdd.volume_fraction[i]
        @. daughter += d * mdd.volume_fraction[i]
    end
    # Calculate ages for each degassing step
    for i in eachindex(parent, daughter)
        age[i] = newton_ar_age(daughter[i], parent[i])
    end
    # Cumulative fraction of parent degassed
    cumsum!(fraction, parent)
    fraction ./= last(fraction)

    return age, fraction
end

function model_ll(mdd::MultipleDomain{T}, σ::T=zero(T)) where {T<:AbstractFloat}
    ll = zero(T)
    for i in eachindex(mdd.age, mdd.age_sigma, mdd.fraction_experimental, mdd.fit)
        if mdd.fit[i]
            model_ageᵢ = linterp1(mdd.model_fraction, mdd.model_age, mdd.fraction_experimental[i])
            ll += norm_ll(mdd.age[i], mdd.age_sigma[i], model_ageᵢ, σ)
        end
    end
    return ll
end

function degassing_ll(mdd::MultipleDomain{T}) where {T<:AbstractFloat}
    ll = zero(T)
    for i in eachindex(mdd.tsteps_experimental, mdd.fraction_experimental, mdd.fit)
        if mdd.fit[i]
            model_fractionᵢ = linterp1(mdd.tsteps_degassing, mdd.model_fraction, mdd.tsteps_experimental[i])
            ll += norm_ll(mdd.fraction_experimental[i], mdd.fraction_experimental_sigma, model_fractionᵢ)
        end
    end
    return ll
end
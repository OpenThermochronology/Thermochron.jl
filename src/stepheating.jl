## -- Functions for modelling experimental degassing schedules

function degas!(mineral::PlanarAr{T}, tsteps_degassing::FloatRange, Tsteps_degassing::AbstractVector{T}, dm::Diffusivity{T}; fuse::Bool=true, redegastracer::Bool=false) where T <: AbstractFloat

    # Constants
    D0 = (dm.D0*10000^2)::T                 # cm^2/sec, converted to micron^2/sec  
    Ea = dm.Ea::T                           # kJ/mol
    R = 0.008314472                         # kJ/(K*mol)
    ΔT = mineral.offset::T + 273.15         # Conversion from C to K, plus temperature offset relative to other samples

    # Calculate effective diffusivity at each time step
    De = mineral.De::Vector{T}
    @assert firstindex(De) == firstindex(Tsteps_degassing)
    @assert lastindex(De) >= lastindex(Tsteps_degassing)
    @turbo for i ∈ eachindex(Tsteps_degassing)
        De[i] = D0 * exp(-Ea / R / (Tsteps_degassing[i] + ΔT)) # micron^2/sec
    end

    # Get time and radius discretization
    dr = step(mineral.rsteps)
    nrsteps = mineral.nrsteps::Int
    ntsteps = length(tsteps_degassing)
    step_tracer = @views(mineral.step_tracer[1:ntsteps])
    step_daughter = @views(mineral.step_daughter[1:ntsteps])
    @assert eachindex(tsteps_degassing) == eachindex(Tsteps_degassing) == Base.OneTo(ntsteps)

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}

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
        # Duration of current timestep
        dt = step_at(tsteps_degassing, i)

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
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No daughter ingrowth on lab timescales
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
    fuse && (step_daughter[ntsteps] = max(daughterᵢ₋, zero(T)))

    # Now diffuse parent isotope tracer (as Ar-39), if neccesary
    if redegastracer || !(0 < sum(step_tracer))

        # Convert diffusivity from that of Ar-40 to that of Ar-39 given Dₗ/Dₕ ~ (mₕ/mₗ)^β 
        # C.f. Luo et al. 2021 (doi: 10.7185/geochemlet.2128) for He in albite melt
        # β=0.355 ± 0.012 at 3000 K, decreasing to β=0.322 ± 0.019 at 1700 K
        De .*= (40/39)^0.3
        
        # Initialize u matrix
        fill!(u, zero(T)) 
        u[2:end-1,1] = mineral.r40K # Model Ar-39 equal to K-40, such that implied ages are already correct
        u[1,1] = u[2,1]
        u[end,1] = 0

        tracerᵢ₋ = nanmean(@views(u[2:end-1, 1]))
        @inbounds for i in Base.OneTo(ntsteps-fuse)
            # Duration of current timestep
            dt = step_at(tsteps_degassing, i)

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
                y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No ingrowth or decay of tracer
            end

            # Invert using tridiagonal matrix algorithm
            # equivalent to u[:,i] = A\y
            lu!(F, A, allowsingular=true)
            ldiv!(F, y)
            u[:,i+1] = y

            tracerᵢ = nanmean(@views(u[2:end-1, i+1]))
            step_tracer[i] = max(tracerᵢ₋ - tracerᵢ, zero(T))
            tracerᵢ₋ = tracerᵢ
        end
        # Degas all remaining tracer in last step if fuse==true
        fuse && (step_tracer[ntsteps] = max(tracerᵢ₋, zero(T)))
    end

    # Return views of the resulting tracer and daughter amounts degassed at each step
    return step_tracer, step_daughter
end
function degas!(mineral::SphericalAr{T}, tsteps_degassing::FloatRange, Tsteps_degassing::AbstractVector{T}, dm::Diffusivity{T}; fuse::Bool=true, redegastracer::Bool=false) where T <: AbstractFloat

    # Constants
    D0 = (dm.D0*10000^2)::T                 # cm^2/sec, converted to micron^2/sec  
    Ea = dm.Ea::T                           # kJ/mol
    R = 0.008314472                         # kJ/(K*mol)
    ΔT = mineral.offset::T + 273.15         # Conversion from C to K, plus temperature offset relative to other samples

    # Calculate effective diffusivity at each time step
    De = mineral.De::Vector{T}
    @assert firstindex(De) == firstindex(Tsteps_degassing)
    @assert lastindex(De) >= lastindex(Tsteps_degassing)
    @turbo for i ∈ eachindex(Tsteps_degassing)
        De[i] = D0 * exp(-Ea / R / (Tsteps_degassing[i] + ΔT)) # micron^2/sec
    end

    # Get time and radius discretization
    rsteps = mineral.rsteps
    dr = step(rsteps)
    nrsteps = mineral.nrsteps::Int
    relvolumes = mineral.relvolumes::Vector{T}
    ntsteps = length(tsteps_degassing)
    step_tracer = @views(mineral.step_tracer[1:ntsteps])
    step_daughter = @views(mineral.step_daughter[1:ntsteps])
    @assert eachindex(tsteps_degassing) == eachindex(Tsteps_degassing) == Base.OneTo(ntsteps)

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = mineral.y

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the Ar profile and r is radius
    u = mineral.u::DenseMatrix{T}
    fill!(u, zero(T)) 
    u[:,1] = y # Start with final Ar profile from end of last geologic inversion

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = mineral.A       # Tridiagonal matrix
    F = mineral.F       # LU object for in-place lu factorization
    
    @inbounds for i in Base.OneTo(ntsteps-fuse)
        # Duration of current timestep
        dt = step_at(tsteps_degassing, i)

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
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No daughter ingrowth on lab timescales
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
    fuse && (step_daughter[ntsteps] = max(daughterᵢ₋, zero(T)))

    # Now diffuse parent isotope tracer (as Ar-39), if neccesary
    if redegastracer || !(0 < sum(step_tracer))
        
        # Convert diffusivity from that of Ar-40 to that of Ar-39 given Dₗ/Dₕ ~ (mₕ/mₗ)^β 
        # C.f. Luo et al. 2021 (doi: 10.7185/geochemlet.2128) for He in albite melt
        # β=0.355 ± 0.012 at 3000 K, decreasing to β=0.322 ± 0.019 at 1700 K
        De .*= (40/39)^0.3

        # Initialize u matrix
        fill!(u, zero(T)) 
        u[2:end-1,1] = mineral.r40K # Model Ar-39 equal to K-40, such that implied ages are already correct
        u[2:end-1,1] .*= rsteps # U-transform for Crank-Nicholson
        u[1,1] = -u[2,1]
        u[end,1] = 0

        @inbounds for i in Base.OneTo(ntsteps-fuse)
            # Duration of current timestep
            dt = step_at(tsteps_degassing, i)

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
                y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No ingrowth or decay of tracer
            end

            # Invert using tridiagonal matrix algorithm
            # equivalent to u[:,i] = A\y
            lu!(F, A, allowsingular=true)
            ldiv!(F, y)
            u[:,i+1] = y
        end

        # Calculate parent (tracer) lost to diffusion at each step
        tracerᵢ₋ = nanmean(@views(u[2:end-1, 1])./=rsteps, relvolumes)
        @inbounds for i in Base.OneTo(ntsteps-fuse)
            tracerᵢ = nanmean(@views(u[2:end-1, i+1])./=rsteps, relvolumes)
            step_tracer[i] = max(tracerᵢ₋ - tracerᵢ, zero(T))
            tracerᵢ₋ = tracerᵢ
        end
        # Degas all remaining tracer in last step if fuse==true
        fuse && (step_tracer[ntsteps] = max(tracerᵢ₋, zero(T)))
    end

    # Return views of the resulting tracer and daughter amounts degassed at each step
    return step_tracer, step_daughter
end
function degas!(mineral::PlanarHe{T}, tsteps_degassing::FloatRange, Tsteps_degassing::AbstractVector{T}, dm::Diffusivity{T}; fuse::Bool=true, redegastracer::Bool=false) where T <: AbstractFloat

    # Constants
    D0 = (dm.D0*10000^2)::T                 # cm^2/sec, converted to micron^2/sec  
    Ea = dm.Ea::T                           # kJ/mol
    R = 0.008314472                         # kJ/(K*mol)
    ΔT = mineral.offset::T + 273.15         # Conversion from C to K, plus temperature offset relative to other samples

    # Calculate effective diffusivity at each time step
    De = mineral.De::Vector{T}
    @assert firstindex(De) == firstindex(Tsteps_degassing)
    @assert lastindex(De) >= lastindex(Tsteps_degassing)
    @turbo for i ∈ eachindex(Tsteps_degassing)
        De[i] = D0 * exp(-Ea / R / (Tsteps_degassing[i] + ΔT)) # micron^2/sec
    end

    # Get time and radius discretization
    dr = step(mineral.rsteps)
    nrsteps = mineral.nrsteps::Int
    ntsteps = length(tsteps_degassing)
    step_tracer = @views(mineral.step_tracer[1:ntsteps])
    step_daughter = @views(mineral.step_daughter[1:ntsteps])
    @assert eachindex(tsteps_degassing) == eachindex(Tsteps_degassing) == Base.OneTo(ntsteps)

    # Common β factor
    β = mineral.β::Vector{T}

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
    
    total_daughter = daughterᵢ₋ = nanmean(@views(u[2:end-1, 1]))
    @inbounds for i in Base.OneTo(ntsteps-fuse)
        # Duration of current timestep
        dt = step_at(tsteps_degassing, i)

        # Update β for current temperature
        # constant across all radii since diffusivity is constant
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
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No daughter ingrowth on lab timescales
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
    fuse && (step_daughter[ntsteps] = max(daughterᵢ₋, zero(T)))

    # Now diffuse isotope tracer (He-3), if neccesary
    if redegastracer || !(0 < sum(step_tracer))

        # Convert diffusivity from that of He-4 to that of He-3 given Dₗ/Dₕ ~ (mₕ/mₗ)^β 
        # C.f. Luo et al. 2021 (doi: 10.7185/geochemlet.2128) for He in albite melt
        # β=0.355 ± 0.012 at 3000 K, decreasing to β=0.322 ± 0.019 at 1700 K
        De .*= (4/3)^0.3
        
        # Initialize u matrix
        fill!(u, zero(T)) 
        u[2:end-1,1] .= total_daughter # Initialize with tracer equal to total daughter, such that results are normalized to Bulk 4He/3He (i.e., Rstep/Rbulk)
        u[1,1] = u[2,1]
        u[end,1] = 0

        tracerᵢ₋ = nanmean(@views(u[2:end-1, 1]))
        @inbounds for i in Base.OneTo(ntsteps-fuse)
            # Duration of current timestep
            dt = step_at(tsteps_degassing, i)

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
                y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No ingrowth or decay of tracer
            end

            # Invert using tridiagonal matrix algorithm
            # equivalent to u[:,i] = A\y
            lu!(F, A, allowsingular=true)
            ldiv!(F, y)
            u[:,i+1] = y

            tracerᵢ = nanmean(@views(u[2:end-1, i+1]))
            step_tracer[i] = max(tracerᵢ₋ - tracerᵢ, zero(T))
            tracerᵢ₋ = tracerᵢ
        end
        # Degas all remaining tracer in last step if fuse==true
        fuse && (step_tracer[ntsteps] = max(tracerᵢ₋, zero(T)))
    end

    # Return views of the resulting tracer and daughter amounts degassed at each step
    return step_tracer, step_daughter
end
function degas!(mineral::SphericalHe{T}, tsteps_degassing::FloatRange, Tsteps_degassing::AbstractVector{T}, dm::Diffusivity{T}; fuse::Bool=true, redegastracer::Bool=false) where T <: AbstractFloat

    # Constants
    D0 = (dm.D0*10000^2)::T                 # cm^2/sec, converted to micron^2/sec  
    Ea = dm.Ea::T                           # kJ/mol
    R = 0.008314472                         # kJ/(K*mol)
    ΔT = mineral.offset::T + 273.15         # Conversion from C to K, plus temperature offset relative to other samples

    # Calculate effective diffusivity at each time step
    De = mineral.De::Vector{T}
    @assert firstindex(De) == firstindex(Tsteps_degassing)
    @assert lastindex(De) >= lastindex(Tsteps_degassing)
    @turbo for i ∈ eachindex(Tsteps_degassing)
        De[i] = D0 * exp(-Ea / R / (Tsteps_degassing[i] + ΔT)) # micron^2/sec
    end

    # Get time and radius discretization
    rsteps = mineral.rsteps
    dr = step(rsteps)
    nrsteps = mineral.nrsteps::Int
    relvolumes = mineral.relvolumes::Vector{T}
    ntsteps = length(tsteps_degassing)
    step_tracer = @views(mineral.step_tracer[1:ntsteps])
    step_daughter = @views(mineral.step_daughter[1:ntsteps])
    @assert eachindex(tsteps_degassing) == eachindex(Tsteps_degassing) == Base.OneTo(ntsteps)

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = mineral.y

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the He profile and r is radius
    u = mineral.u::DenseMatrix{T}
    fill!(u, zero(T)) 
    u[:,1] = y # Start with final (coordinate transform'd) He profile from end of last geologic inversion

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = mineral.A       # Tridiagonal matrix
    F = mineral.F       # LU object for in-place lu factorization
    
    @inbounds for i in Base.OneTo(ntsteps-fuse)
        # Duration of current timestep
        dt = step_at(tsteps_degassing, i)

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
            y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No daughter ingrowth on lab timescales
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i+1] = y
    end

    # Calculate daughter lost to diffusion at each step
    total_daughter = daughterᵢ₋ = nanmean(@views(u[2:end-1, 1])./=rsteps, relvolumes)
    @inbounds for i in Base.OneTo(ntsteps-fuse)
        daughterᵢ = nanmean(@views(u[2:end-1, i+1])./=rsteps, relvolumes)
        step_daughter[i] = max(daughterᵢ₋ - daughterᵢ, zero(T))
        daughterᵢ₋ = daughterᵢ
    end
    # Degas all remaining daughter in last step if fuse==true
    fuse && (step_daughter[ntsteps] = max(daughterᵢ₋, zero(T)))

    # Now diffuse isotope tracer (He-3), if neccesary
    if redegastracer || !(0 < sum(step_tracer))
        
        # Convert diffusivity from that of He-4 to that of He-3 given Dₗ/Dₕ ~ (mₕ/mₗ)^β 
        # C.f. Luo et al. 2021 (doi: 10.7185/geochemlet.2128) for He in albite melt
        # β=0.355 ± 0.012 at 3000 K, decreasing to β=0.322 ± 0.019 at 1700 K
        De .*= (4/3)^0.3

        # Initialize u matrix
        fill!(u, zero(T)) 
        u[2:end-1,1] .= total_daughter # Initialize with tracer equal to total daughter, such that results are normalized to Bulk 4He/3He (i.e., Rstep/Rbulk)
        u[2:end-1,1] .*= rsteps # U-transform for Crank-Nicholson
        u[1,1] = -u[2,1]
        u[end,1] = 0

        @inbounds for i in Base.OneTo(ntsteps-fuse)
            # Duration of current timestep
            dt = step_at(tsteps_degassing, i)

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
                y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No ingrowth or decay of tracer
            end

            # Invert using tridiagonal matrix algorithm
            # equivalent to u[:,i] = A\y
            lu!(F, A, allowsingular=true)
            ldiv!(F, y)
            u[:,i+1] = y
        end

        # Calculate tracer lost to diffusion at each step
        tracerᵢ₋ = nanmean(@views(u[2:end-1, 1])./=rsteps, relvolumes)
        @inbounds for i in Base.OneTo(ntsteps-fuse)
            tracerᵢ = nanmean(@views(u[2:end-1, i+1])./=rsteps, relvolumes)
            step_tracer[i] = max(tracerᵢ₋ - tracerᵢ, zero(T))
            tracerᵢ₋ = tracerᵢ
        end
        # Degas all remaining tracer in last step if fuse==true
        fuse && (step_tracer[ntsteps] = max(tracerᵢ₋, zero(T)))
    end

    # Return views of the resulting tracer and daughter amounts degassed at each step
    return step_tracer, step_daughter
end

## --- Age and likelihood functions for step heating data

function modelage(mdd::MultipleDomain{T}, Tsteps::AbstractVector, dm::MDDiffusivity{T}; redegastracer::Bool=false) where {T<:AbstractFloat}
    age = fill!(mdd.model_age, zero(T))
    tracer = fill!(mdd.model_tracer, zero(T))
    daughter = fill!(mdd.model_daughter, zero(T))
    fraction = fill!(mdd.model_fraction, zero(T))
    fuse = mdd.fuse::Bool
    # Degas
    for i in eachindex(mdd.domains, mdd.volume_fraction)
        domain = mdd.domains[i]
        modelage(domain, Tsteps, dm[i])
        p, d = degas!(domain, mdd.tsteps_degassing, mdd.Tsteps_degassing, dm[i]; fuse, redegastracer)
        @. tracer += p * mdd.volume_fraction[i]
        @. daughter += d * mdd.volume_fraction[i]
    end
    # Calculate ages for each degassing step
    for i in eachindex(tracer, daughter)
        age[i] = newton_ar_age(daughter[i], tracer[i])
    end
    # Cumulative fraction of tracer degassed
    cumsum!(fraction, tracer)
    fraction ./= last(fraction)

    return age, fraction
end

function model_ll(mdd::MultipleDomain{T}, σ::T=zero(T); rescalemdd=true) where {T<:AbstractFloat}
    ll = zero(T)
    for i in eachindex(mdd.age, mdd.age_sigma, mdd.midpoint_experimental, mdd.fit)
        if mdd.fit[i]
            model_ageᵢ = linterp1(mdd.model_fraction, mdd.model_age, mdd.midpoint_experimental[i])
            ll += norm_ll(mdd.age[i], mdd.age_sigma[i], model_ageᵢ, σ)
        end
    end
    rescalemdd && (ll /= sqrt(count(mdd.fit)))
    return ll
end

function degassing_ll(mdd::MultipleDomain{T}; rescalemdd=true) where {T<:AbstractFloat}
    ll = zero(T)
    fit_until = findlast(mdd.fit)
    for i in eachindex(mdd.tsteps_experimental, mdd.fraction_experimental, mdd.fit)
        if i <= fit_until
            model_fractionᵢ = linterp1(mdd.tsteps_degassing, mdd.model_fraction, mdd.tsteps_experimental[i])
            ll += norm_ll(mdd.fraction_experimental[i], mdd.fraction_experimental_sigma, model_fractionᵢ)
        end
    end
    rescalemdd && (ll /= sqrt(fit_until))
    return ll
end
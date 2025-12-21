## --- Calculate β parameters in-place with current diffusivity information

function updatebeta!(β::Vector{T}, mineral::ZirconHe{T}, dm::ZRDAAM{T}, dt::T, TK::T, tstep::Int, diffusivityratio=one(T); setting::Symbol=:geological) where {T}
    # The annealed damage matrix
    # We will use the last timestep of the geological annealed damage matrix,
    # Assuming He diffuses out faster than annealing occurs during laboratory degassing
    annealeddamage = mineral.annealeddamage::Matrix{T}

    # Damage and annealing constants
    DzEa = dm.DzEa::T                           # [kJ/mol]
    DN17Ea = dm.DN17Ea::T                       # [kJ/mol]
    lint0 = dm.lint0::T                         # [nm]
    SV = dm.SV::T                               # [1/nm]
    Bα = dm.Bα::T                               # [g/alpha] mass of amorphous material produced per alpha decay
    Phi = dm.Phi::T                             # [unitless]
    R = 0.008314472                             # [kJ/(K*mol)]
    if setting === :laboratory
        DzD0 = (dm.DzD0*10000^2)::T             # [micron^2/sec], converted from [cm^2/sec]
        DN17D0 = (dm.DN17D0*10000^2)::T         # [micron^2/sec], converted from [cm^2/sec]
        damagestep = lastindex(annealeddamage, 1)
    else
        DzD0 = (dm.DzD0*10000^2*SEC_MYR)::T     # [micron^2/Myr], converted from [cm^2/sec]
        DN17D0 = (dm.DN17D0*10000^2*SEC_MYR)::T # [micron^2/Myr], converted from [cm^2/sec]
        damagestep = tstep
    end

    # Radial step size
    dr = step(mineral.rsteps)::T
    nrsteps = mineral.nrsteps
    @assert eachindex(β) == Base.OneTo(nrsteps)

    # Endmember diffusion constants
    Dz = DzD0 * exp(-DzEa / (R * TK)) * diffusivityratio # [micron^2/sec]
    DN17 = DN17D0 * exp(-DN17Ea / (R * TK)) * diffusivityratio # [micron^2/sec

    # Each radial step except first and latst
    @turbo for k in Base.OneTo(nrsteps-2)
        dam = annealeddamage[damagestep,k]
        fₐ = 1-exp(-Bα*dam*Phi)
        τ = (lint0/(4.2 / ((1-exp(-Bα*dam)) * SV) - 2.5))^2
        De = 1 / ((1-fₐ)^3 / (Dz/τ) + fₐ^3 / DN17)  # Effetive diffusivity
        β[k+1] = 2 * dr^2 / (De*dt) # Shifted by 1 because β[1] is implicit point at negative radius
    end
    # First and last radial step
    β[1] = β[2]
    β[end] = β[end-1]

    return β
end
function updatebeta!(β::Vector{T}, mineral::ApatiteHe{T}, dm::RDAAM{T}, dt::T, TK::T, tstep::Int, diffusivityratio=one(T); setting::Symbol=:laboratory) where {T}
    # The annealed damage matrix
    # We will use the last timestep of the geological annealed damage matrix,
    # Assuming He diffuses out faster than annealing occurs during laboratory degassing
    annealeddamage = mineral.annealeddamage::Matrix{T}

    # Damage and annealing constants
    EaL = dm.EaL::T                         # [kJ/mol]
    EaTrap = dm.EaTrap::T                   # [kJ/mol]
    etaq = dm.etaq::T                       # Durango ηq
    psi = dm.psi::T                         # [unitless]
    omega = dm.omega::T                     # [unitless]
    rhoap = dm.rhoap::T                     # [g/cm^3]
    L = dm.L::T                             # [cm]
    lambdaf = dm.lambdaf::T                 # [1/time]
    lambdaD = dm.lambdaD::T                 # [1/time]
    R = 0.008314472                         # [kJ/(K*mol)]
    if setting === :laboratory
        D0L = (dm.D0L*10000^2)::T           # [micron^2/sec], converted from [cm^2/sec]  
        damagestep = lastindex(annealeddamage, 1)
    else
        D0L = (dm.D0L*10000^2*SEC_MYR)::T   # [micron^2/Myr], converted from [cm^2/sec]  
        damagestep = tstep
    end

    # Radial step size
    dr = step(mineral.rsteps)::T
    nrsteps = mineral.nrsteps
    @assert eachindex(β) == Base.OneTo(nrsteps)

    # Conversion factor from alphas/g to track length cm/cm^3
    damage_conversion = rhoap*(lambdaf/lambdaD)*etaq*L

    # Normal and trapping diffusivities
    DL = D0L * exp(-EaL / (R * TK)) * diffusivityratio # [micron^2/t]
    Dtrap = exp( EaTrap / (R * TK)) # [unitless]

    # Each radial step except first and latst
    @turbo for k in Base.OneTo(nrsteps-2)
        track_density = annealeddamage[damagestep, k]*damage_conversion # [cm/cm3]
        trap = (psi*track_density + omega*track_density^3)*Dtrap
        De = DL/(trap+1) # [micron^2/t]
        β[k+1] = 2 * dr^2 / (De*dt)
    end
    # First and last radial step
    β[1] = β[2]
    β[end] = β[end-1]

    return β
end
function updatebeta!(β::Vector{T}, mineral::NobleGasSample{T}, dm::Diffusivity{T}, dt::T, TK::T, ::Int, diffusivityratio=one(T); setting::Symbol=:laboratory) where {T}
    # Constants
    Ea = dm.Ea::T                           # [kJ/mol]
    R = 0.008314472                         # [kJ/(K*mol)]
    if setting === :laboratory
        D0 = (dm.D0*10000^2)::T             # [micron^2/sec], converted from [cm^2/sec]
    else
        D0 = (dm.D0*10000^2*SEC_MYR)::T     # [micron^2/Myr], converted from [cm^2/sec]
    end

    # Radial step size
    dr = step(mineral.rsteps)::T

    # Diffusivity
    De = D0 * exp(-Ea / (R * TK)) * diffusivityratio # [micron^2/t]

    # Update β for current temperature
    fill!(β, 2 * dr^2 / (De*dt))

    return β
end

# --- Core Crank-Nicolson solvers

function crank_nicolson!(mineral::PlanarNobleGas{T}, tsteps::AbstractVector{T}, Tsteps::AbstractVector{T}, dm::Diffusivity{T}; 
        fuse::Bool=false, 
        partitiondaughter::Bool=false, 
        diffusivityratio::Number=one(T), 
        setting::Symbol=:laboratory,
    ) where {T <: AbstractFloat}

    @assert (setting===:laboratory) || (setting===:geological)
    partitiondaughter &= (setting===:geological)
    ingrowth = (setting===:geological)

    # Temperature conversion
    ΔT = T(273.15)
    (setting === :geological) && (ΔT += mineral.offset)

    # Check time and radius discretization
    nrsteps = mineral.nrsteps::Int
    ntsteps = length(tsteps)
    @assert eachindex(tsteps) == eachindex(Tsteps) == Base.OneTo(ntsteps)

    # Variables related to daughter ingrowth/deposition
    bulkdaughter = zero(T)
    bulkdeposition = mineral.bulkdeposition::Vector{T}
    deposition = mineral.deposition::Matrix{T}
    @assert eachindex(bulkdeposition) == axes(deposition, 1) == Base.OneTo(ntsteps)

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the Ar profile and r is radius
    u = mineral.u::DenseMatrix{T}
    @assert axes(u,2) == 1:ntsteps+1
    @assert axes(u,1) == 1:nrsteps

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = mineral.y

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = mineral.A       # Tridiagonal matrix
    F = mineral.F       # LU object for in-place lu factorization

    @inbounds for i in Base.OneTo(ntsteps-fuse)
        # Duration and temperature of current timestep
        dt = step_at(tsteps, i)
        TK = Tsteps[i]+ΔT

        # Beta factor containing diffusivity information for each radial step
        updatebeta!(β, mineral, dm, dt, TK, i, diffusivityratio; setting)

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

        # Optionally set external boundary condition assuming daughter partitioning between
        # grain and intragranular medium, following the partition coefficient measurements
        # of Baxter, Asimow, and Farley (2006) https://doi.org/10.1016/j.gca.2006.09.011
        if partitiondaughter
            bulkdaughter += bulkdeposition[i]
            y[nrsteps] = bulkdaughter * fraction_internal(TK, mineral)
        end

        # RHS of tridiagonal Crank-Nicolson equation for regular grid cells.
        # From Ketcham (2005) https://doi.org/10.2138/rmg.2005.58.11
        if ingrowth
            @turbo for k = 2:nrsteps-1
                𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
                y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊ - deposition[i, k-1]*β[k]
            end
        else
            @turbo for k = 2:nrsteps-1
                𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
                y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No diffusant depositon on lab timescales
            end
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i+1] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i+1] = y
    end

    return mineral
end
function crank_nicolson!(mineral::SphericalNobleGas{T}, tsteps::AbstractVector{T}, Tsteps::AbstractVector{T}, dm::DiffusivityModel{T}; 
        fuse::Bool=false, 
        partitiondaughter::Bool=false, 
        diffusivityratio::Number=one(T), 
        setting::Symbol=:laboratory,
    ) where {T <: AbstractFloat}

    @assert (setting===:laboratory) || (setting===:geological)
    partitiondaughter &= (setting===:geological)
    ingrowth = (setting===:geological)

    # Temperature conversion
    ΔT = T(273.15)
    (setting === :geological) && (ΔT += mineral.offset)

    # Check time and radius discretization
    nrsteps = mineral.nrsteps::Int
    rsteps = mineral.rsteps::AbstractVector{T}
    @assert eachindex(rsteps) == Base.OneTo(nrsteps-2)
    ntsteps = length(tsteps)
    @assert eachindex(tsteps) == eachindex(Tsteps) == Base.OneTo(ntsteps)

    # Variables related to daughter ingrowth/deposition
    bulkdaughter = zero(T)
    bulkradius = last(rsteps) + step(rsteps)
    bulkdeposition = mineral.bulkdeposition::Vector{T}
    deposition = mineral.deposition::Matrix{T}
    @assert eachindex(bulkdeposition) == axes(deposition, 1) == Base.OneTo(ntsteps)

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the Ar profile and r is radius
    u = mineral.u::DenseMatrix{T}
    @assert axes(u,2) == 1:ntsteps+1
    @assert axes(u,1) == 1:nrsteps

    # Common β factor is constant across all radii since diffusivity is constant
    β = mineral.β::Vector{T}

    # Vector for RHS of Crank-Nicholson equation with regular grid cells
    y = mineral.y

    # Tridiagonal matrix for LHS of Crank-Nicolson equation with regular grid cells
    A = mineral.A       # Tridiagonal matrix
    F = mineral.F       # LU object for in-place lu factorization

    @inbounds for i in Base.OneTo(ntsteps-fuse)
        # Duration and temperature of current timestep
        dt = step_at(tsteps, i)
        TK = Tsteps[i]+ΔT

        # Beta factor containing diffusivity information for each radial step
        updatebeta!(β, mineral, dm, dt, TK, i, diffusivityratio; setting)

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

        # Optionally set external boundary condition assuming daughter partitioning between
        # grain and intragranular medium, following the partition coefficient measurements
        # of Baxter, Asimow, and Farley (2006) https://doi.org/10.1016/j.gca.2006.09.011
        if partitiondaughter
            bulkdaughter += bulkdeposition[i]
            y[nrsteps] = bulkradius * bulkdaughter * fraction_internal(TK, mineral)
        end

        # RHS of tridiagonal Crank-Nicolson equation for regular grid cells.
        # From Ketcham (2005) https://doi.org/10.2138/rmg.2005.58.11
        if ingrowth
            @turbo for k = 2:nrsteps-1
                𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
                y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊ - rsteps[k-1]*deposition[i, k-1]*β[k]
            end
        else
            @turbo for k = 2:nrsteps-1
                𝑢ⱼ, 𝑢ⱼ₋, 𝑢ⱼ₊ = u[k, i], u[k-1, i], u[k+1, i]
                y[k] = (2.0-β[k])*𝑢ⱼ - 𝑢ⱼ₋ - 𝑢ⱼ₊    # No diffusant depositon on lab timescales
            end
        end

        # Invert using tridiagonal matrix algorithm
        # equivalent to u[:,i+1] = A\y
        lu!(F, A, allowsingular=true)
        ldiv!(F, y)
        u[:,i+1] = y
    end

    return mineral
end

## --- How much diffusant was lost during a Crank-Nicolson run?

function diffusant_lost!(step_diffusant::AbstractVector{T}, mineral::PlanarNobleGas{T}; fuse::Bool=true) where {T}
    # Check time and radius discretization
    ntsteps = length(step_diffusant)
    @assert eachindex(step_diffusant) == Base.OneTo(ntsteps)

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the Ar profile and r is radius
    u = mineral.u::DenseMatrix{T}
    @assert axes(u,2) == Base.OneTo(ntsteps+1)
    @assert axes(u,1) == Base.OneTo(mineral.nrsteps)

    # Calculate diffusant lost to diffusion at each step
    total_diffusant = diffusantᵢ₋ = nanmean(@views(u[2:end-1, 1]))
    @inbounds for i in Base.OneTo(ntsteps-fuse)
        diffusantᵢ = nanmean(@views(u[2:end-1, i+1]))
        step_diffusant[i] = max(diffusantᵢ₋ - diffusantᵢ, zero(T))
        diffusantᵢ₋ = diffusantᵢ
    end
    # Degas all remaining diffusant in last step if fuse==true
    fuse && (step_diffusant[ntsteps] = max(diffusantᵢ₋, zero(T)))

    return total_diffusant::T
end
function diffusant_lost!(step_diffusant::AbstractVector{T}, mineral::SphericalNobleGas{T}; fuse::Bool=true) where {T}
    # Check time and radius discretization
    ntsteps = length(step_diffusant)
    @assert eachindex(step_diffusant) == Base.OneTo(ntsteps)
    nrsteps = mineral.nrsteps::Int
    rsteps = mineral.rsteps::AbstractVector{T}
    relvolumes = mineral.relvolumes::Vector{T}
    @assert eachindex(relvolumes) == eachindex(rsteps) == Base.OneTo(nrsteps-2)

    # Output matrix for all timesteps
    # u = v*r is the coordinate transform (u-substitution) for the Crank-
    # Nicholson equations where v is the Ar profile and r is radius
    u = mineral.u::DenseMatrix{T}
    @assert axes(u,2) == Base.OneTo(ntsteps+1)
    @assert axes(u,1) == Base.OneTo(nrsteps)

    # Calculate diffusant lost to diffusion at each step
    total_diffusant = diffusantᵢ₋ = nanmean(@views(u[2:end-1, 1])./=rsteps, relvolumes)
    @inbounds for i in Base.OneTo(ntsteps-fuse)
        diffusantᵢ = nanmean(@views(u[2:end-1, i+1])./=rsteps, relvolumes)
        step_diffusant[i] = max(diffusantᵢ₋ - diffusantᵢ, zero(T))
        diffusantᵢ₋ = diffusantᵢ
    end
    # Degas all remaining diffusant in last step if fuse==true
    fuse && (step_diffusant[ntsteps] = max(diffusantᵢ₋, zero(T)))

    return total_diffusant::T
end

## --- How much diffusant is left at the end of a Crank-Nicolson run?

function final_diffusant(mineral::PlanarNobleGas{T}) where {T<:AbstractFloat}
    vfinal = @views mineral.u[2:end-1,end]
    return nanmean(vfinal)::T # Atoms/gram
end
function final_diffusant(mineral::SphericalNobleGas{T}) where {T<:AbstractFloat}
    # Convert from u (coordinate-transform'd conc.) to v (conc.)
    vfinal = @views mineral.u[2:end-1,end]
    vfinal ./= mineral.rsteps
    return nanmean(vfinal, mineral.relvolumes)::T # Atoms/gram
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

## --- End of File

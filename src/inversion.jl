
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
    function maxdiff(x::AbstractVector)
        i₀ = firstindex(x)
        δₘ = x[i₀+1] - x[i₀]
        if length(x) > 2
            last = x[i₀+1]
            @inbounds for i ∈ (i₀+2):(i₀+length(x)-1)
                δᵢ = x[i] - last
                if δᵢ > δₘ
                    δₘ = δᵢ
                end
                last = x[i]
            end
        end
        return δₘ
    end
    function mindiff(x::AbstractVector)
        i₀ = firstindex(x)
        δₘ = x[i₀+1] - x[i₀]
        if length(x) > 2
            last = x[i₀+1]
            @inbounds for i ∈ (i₀+2):(i₀+length(x)-1)
                δᵢ = x[i] - last
                if δᵢ < δₘ
                    δₘ = δᵢ
                end
                last = x[i]
            end
        end
        return δₘ
    end
    function maxabsdiff(x::AbstractVector)
        i₀ = firstindex(x)
        δₘ = abs(x[i₀+1] - x[i₀])
        if length(x) > 2
            last = x[i₀+1]
            @inbounds for i ∈ (i₀+2):(i₀+length(x)-1)
                δᵢ = abs(x[i] - last)
                if δᵢ > δₘ
                    δₘ = δᵢ
                end
                last = x[i]
            end
        end
        return δₘ
    end

    # Check if point k is distinct from other points in list within ± δ
    function isdistinct(points::DenseArray, npoints::Int, k::Int, δ::Number)
        @inbounds for i = 1:npoints
            if i!=k && abs(points[i] - points[k]) < δ
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

    """
    ```julia
    MCMC_vartcryst(data::NamedTuple, model::NamedTuple, npoints::Int, agepoints::Vector, Tpoints::Vector, unconf::Unconformity, boundary::Boundary)
    ```
    Markov chain Monte Carlo time-Temperature inversion of the data specified in `data` and model parameters specified by `model`.

    Returns a tuple of distributions `(tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist)`

    ## Examples
    ```julia
    tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist = MCMC_vartcryst(data, model, npoints, agepoints, Tpoints, unconf, boundary)
    ```
    """
    function MCMC_vartcryst(data::NamedTuple, model::NamedTuple, npoints::Int, agepoints::DenseVector{T}, Tpoints::DenseVector{T}, boundary::Boundary{T}, unconf::Unconformity{T}) where T <: AbstractFloat
        # Sanitize inputs
        @assert firstindex(agepoints) === 1
        @assert firstindex(Tpoints) === 1
        halfwidth = T.(data.halfwidth)::DenseVector{T}
        halfwidth = T.(data.halfwidth)::DenseVector{T}
        crystAge = T.(data.crystAge)::DenseVector{T}
        HeAge = T.(data.HeAge)::DenseVector{T}
        HeAge_sigma = T.(data.HeAge_sigma)::DenseVector{T}
        U = T.(data.U)::DenseVector{T}
        Th = T.(data.Th)::DenseVector{T}
        nsteps = model.nsteps::Int
        maxpoints = model.maxpoints::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        totalpoints = maxpoints + boundary.npoints + unconf.npoints::Int
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(model.dTmax)::T
        agesteps = T.(model.agesteps)::DenseVector{T}
        tsteps = T.(model.tsteps)::DenseVector{T}
        tinit = T(model.tinit)::T
        Tinit = T(model.Tinit)::T
        Tnow = T(model.Tnow)::T
        dt = T(model.dt)::T
        dr = T(model.dr)::T
        σmodel = T(model.σmodel)::T
        σannealing = T(model.σannealing)::T
        λannealing = T(model.λannealing)::T

        # Calculate number of boundary and unconformity points and allocate buffer for interpolating
        agepointbuffer = similar(agepoints, totalpoints)::DenseVector{T}
        Tpointbuffer = similar(agepoints, totalpoints)::DenseVector{T}
        knot_index = similar(agesteps, Int)::DenseVector{Int}

        # Calculate model ages for initial proposal
        ages = collectto!(agepointbuffer, view(agepoints, 1:npoints), boundary.agepoints, unconf.agepoints)::StridedVector{T}
        temperatures = collectto!(Tpointbuffer, view(Tpoints, 1:npoints), boundary.Tpoints, unconf.Tpoints)::StridedVector{T}
        Tsteps = linterp1s(ages, temperatures, agesteps)::DenseVector{T}
        calcHeAges = similar(HeAge)::DenseVector{T}
        pr, Teq = anneal(dt, tsteps, Tsteps, ZRDAAM()) # Damage annealing history
        pr::DenseMatrix{T}
        Teq::DenseVector{T}

        # Prepare a Mineral object for each analysis
        diffusionmodel = (DzEa=T(model.DzEa), DzD0=T(model.DzD0), DN17Ea=T(model.DN17Ea), DN17D0=T(model.DN17D0))
        zircons = Array{Zircon{T}}(undef, length(halfwidth))::Vector{Zircon{T}}
        for i=1:length(zircons)
            # Iterate through each grain, calculate the modeled age for each
            first_index = 1 + floor(Int64,(tinit - crystAge[i])/dt)
            zircons[i] = Zircon(halfwidth[i], dr, U[i], Th[i], dt, agesteps[first_index:end])
            calcHeAges[i] = HeAgeSpherical(zircons[i], @views(Tsteps[first_index:end]), @views(pr[first_index:end,first_index:end]), diffusionmodel)::T
        end

        # Simulated annealing of uncertainty
        σₐ = simannealsigma.(1, HeAge_sigma, σmodel, σannealing, λannealing)::DenseVector{T}
        σ = sqrt.(HeAge_sigma.^2 .+ σmodel^2)

        # Log-likelihood for initial proposal
        ll = normpdf_ll(HeAge, σₐ, calcHeAges)
        if simplified
            ll -= log(npoints)
        end

        # Variables to hold proposals
        llₚ = ll
        npointsₚ = npoints
        agepointsₚ = copy(agepoints)::DenseVector{T}
        Tpointsₚ = copy(Tpoints)::DenseVector{T}
        calcHeAgesₚ = copy(calcHeAges)::DenseVector{T}
        Tstepsₚ = copy(Tsteps)::DenseVector{T}

        # distributions to populate
        tpointdist = zeros(T, totalpoints, nsteps)
        Tpointdist = zeros(T, totalpoints, nsteps)
        HeAgedist = zeros(T, length(HeAge), nsteps)
        σⱼtdist = zeros(T, nsteps)
        σⱼTdist = zeros(T, nsteps)
        lldist = zeros(T, nsteps)
        ndist = zeros(Int, nsteps)
        acceptancedist = zeros(Bool, nsteps)

        # Standard deviations of Gaussian proposal ("jumping") distributions
        # for temperature and time
        σⱼt = tinit/60
        σⱼT = Tinit/60
        k = 1 # Index of chosen t-T point

        # Proposal probabilities (must sum to 1)
        move = 0.64
        birth = 0.15 # Must equal death
        death = 0.15 # Must equal birth
        movebounds = 0.06

        # Number of times to attempt to satisfy reheating rate: # Should be large
        # enough that proposal probabilities are unchanged, but low enough to prevent
        # an infinite loop
        nattempts = 100_000

        progress = Progress(nsteps, dt=1, desc="Running MCMC ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps

            # Copy proposal from last accepted solution
            npointsₚ = npoints
            copyto!(agepointsₚ, agepoints)
            copyto!(Tpointsₚ, Tpoints)
            copyto!(unconf.agepointsₚ, unconf.agepoints)
            copyto!(unconf.Tpointsₚ, unconf.Tpoints)
            copyto!(boundary.Tpointsₚ, boundary.Tpoints)

            # Adjust the proposal
            r = rand()
            if r < move
                # Move one t-T point
                k = ceil(Int, rand() * npoints)
                for attempt ∈ 1:nattempts

                # Move the age of one model point
                agepointsₚ[k] += randn() * σⱼt
                # if agepointsₚ[k] < (0 + dt)
                #     # Reflecting boundary condition at (0 + dt)
                #     agepointsₚ[k] = (0 + dt) - (agepointsₚ[k] - (0 + dt))
                # elseif agepointsₚ[k] > (tinit - dt)
                #     # Reflecting boundary condition at (tinit - dt)
                #     agepointsₚ[k] = (tinit - dt) - (agepointsₚ[k] - (tinit - dt))
                # end

                # Move the Temperature of one model point
                Tpointsₚ[k] += randn() * σⱼT
                # if Tpointsₚ[k] < Tnow
                #     # Reflecting boundary conditions at Tnow (0)
                #     Tpointsₚ[k] = Tnow - (Tpointsₚ[k] - Tnow)
                # elseif Tpointsₚ[k] > Tinit
                #     # Reflecting boundary conditions at Tinit
                #     Tpointsₚ[k] = Tinit - (Tpointsₚ[k] - Tinit)
                # end

                # Circular boundary conditions
                agepointsₚ[k] = mod(agepointsₚ[k]-dt, tinit-2dt) + dt
                Tpointsₚ[k] = mod(Tpointsₚ[k]-Tnow, Tinit-Tnow) + Tnow
                # agepointsₚ[k] = min(max(agepointsₚ[k], dt), tinit)
                # Tpointsₚ[k] = min(max(Tpointsₚ[k], Tnow), Tinit)

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, unconf.agepointsₚ)::StridedVector{T}
                temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, unconf.Tpointsₚ)::StridedVector{T}
                linterp1s!(Tstepsₚ, knot_index, ages, temperatures, agesteps)

                # Retry unless we have satisfied the maximum reheating rate
                if isdistinct(agepointsₚ, npointsₚ, k, 2dt) && maxdiff(Tstepsₚ) < dTmax
                    break
                end
                if attempt == nattempts
                    @warn """`move` proposals failed to satisfy reheating rate limit
                    maxdiff: $(maxdiff(Tstepsₚ))
                    ages: $(agepointsₚ[1:npointsₚ])
                    temperatures: $(Tpointsₚ[1:npointsₚ])
                    σⱼt: $(σⱼt)
                    σⱼT: $(σⱼT)"""
                end
                # Copy last accepted solution to re-modify if we don't break
                copyto!(agepointsₚ, agepoints)
                copyto!(Tpointsₚ, Tpoints)
                end
            elseif (r < move+birth) && (npointsₚ < maxpoints)
                # Birth: add a new model point
                npointsₚ = npoints + 1
                for attempt ∈ 1:nattempts
                agepointsₚ[npointsₚ] = 0 + rand()*(tinit-0)
                Tpointsₚ[npointsₚ] = Tnow + rand()*(Tinit-Tnow)

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, unconf.agepointsₚ)::StridedVector{T}
                temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, unconf.Tpointsₚ)::StridedVector{T}
                linterp1s!(Tstepsₚ, knot_index, ages, temperatures, agesteps)

                # Retry unless we have satisfied the maximum reheating rate
                if isdistinct(agepointsₚ, npointsₚ, npointsₚ, 2dt) && maxdiff(Tstepsₚ) < dTmax
                    break
                end
                if attempt == nattempts
                    @warn """new point proposals failed to satisfy reheating rate limit
                    maxdiff: $(maxdiff(Tstepsₚ))
                    ages: $(agepointsₚ[1:npointsₚ])
                    temperatures: $(Tpointsₚ[1:npointsₚ])
                    σⱼt: $(σⱼt)
                    σⱼT: $(σⱼT)"""
                end
                end
            elseif (r < move+birth+death) && (r >= move+birth) && (npointsₚ > minpoints)
                # Death: remove a model point
                npointsₚ = npoints - 1 # Delete last point in array from proposal
                for attempt ∈ 1:nattempts
                k = ceil(Int, rand()*npoints) # Choose point to delete
                agepointsₚ[k] = agepointsₚ[npoints]
                Tpointsₚ[k] = Tpointsₚ[npoints]

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, unconf.agepointsₚ)::StridedVector{T}
                temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, unconf.Tpointsₚ)::StridedVector{T}
                linterp1s!(Tstepsₚ, knot_index, ages, temperatures, agesteps)

                # Retry unless we have satisfied the maximum reheating rate
                (maxdiff(Tstepsₚ) < dTmax) && break
                if attempt == nattempts
                    @warn """point removal proposals failed to satisfy reheating rate limit
                    maxdiff: $(maxdiff(Tstepsₚ))
                    ages: $(agepointsₚ[1:npointsₚ])
                    temperatures: $(Tpointsₚ[1:npointsₚ])
                    σⱼt: $(σⱼt)
                    σⱼT: $(σⱼT)"""
                end
                # Copy last accepted solution to re-modify if we don't break
                copyto!(agepointsₚ, agepoints)
                copyto!(Tpointsₚ, Tpoints)
                end
            else
                # Move boundary conditions
                for attempt ∈ 1:nattempts
                # Move the temperatures of the starting or ending boundaries
                # @. boundary.Tpointsₚ = boundary.T₀ + rand()*boundary.ΔT
                k = rand(1:boundary.npoints)
                boundary.Tpointsₚ[k] = boundary.T₀[k] + rand()*boundary.ΔT[k]

                # If there's an imposed unconformity, adjust within parameters
                if unconf.npoints > 0
                    @. unconf.agepointsₚ = unconf.Age₀ + rand()*unconf.ΔAge
                    @. unconf.Tpointsₚ = unconf.T₀ + rand()*unconf.ΔT
                end

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, unconf.agepointsₚ)::StridedVector{T}
                temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, unconf.Tpointsₚ)::StridedVector{T}
                linterp1s!(Tstepsₚ, knot_index, ages, temperatures, agesteps)

                # Retry unless we have satisfied the maximum reheating rate
                (maxdiff(Tstepsₚ) < dTmax) && break
                if attempt == nattempts
                    @warn """`movebounds` proposals failed to satisfy reheating rate limit
                    maxdiff: $(maxdiff(Tstepsₚ))
                    ages: $(agepointsₚ[1:npointsₚ])
                    temperatures: $(Tpointsₚ[1:npointsₚ])
                    σⱼt: $(σⱼt)
                    σⱼT: $(σⱼT)"""
                end
                # Copy last accepted solution to re-modify if we don't break
                copyto!(unconf.agepointsₚ, unconf.agepoints)
                copyto!(unconf.Tpointsₚ, unconf.Tpoints)
                copyto!(boundary.Tpointsₚ, boundary.Tpoints)
                end
            end

            if any(isnan, view(agepointsₚ, 1:npointsₚ)) ||  any(isnan, view(Tpointsₚ, 1:npointsₚ))
                @warn """`NaN`s detected!
                ages: $(agepointsₚ[1:npointsₚ])
                temperatures: $(Tpointsₚ[1:npointsₚ])
                σⱼt: $(σⱼt)
                σⱼT: $(σⱼT)
                r: $r
                """
            end

             # Calculate model ages for each grain
            anneal!(pr, Teq, dt, tsteps, Tstepsₚ, ZRDAAM())
            for i=1:length(zircons)
                first_index = 1 + floor(Int64,(tinit - crystAge[i])/dt)
                calcHeAgesₚ[i] = HeAgeSpherical(zircons[i], @views(Tstepsₚ[first_index:end]), @views(pr[first_index:end,first_index:end]), diffusionmodel)::T
            end

            # Calculate log likelihood of proposal
            σₐ .= simannealsigma.(n, HeAge_sigma, σmodel, σannealing, λannealing)
            llₚ = normpdf_ll(HeAge, σₐ, calcHeAgesₚ)
            llₗ = normpdf_ll(HeAge, σₐ, calcHeAges) # Recalulate last one too with new σₐ
            if simplified # slightly penalize more complex t-T paths
                llₚ -= log(npointsₚ)
                llₗ -= log(npoints)
            end

            # Accept or reject proposal based on likelihood
            # To avoid numerical problems with diffusion code, also reject proposal
            # if maximum proposed heating rate is greater than 25C per timestep.
            # (Fast cooling should not be a problem, however)
            if log(rand()) < (llₚ - llₗ)

                # Update jumping distribution based on size of current accepted move
                if dynamicjumping && r < move
                    if agepointsₚ[k] != agepoints[k]
                        σⱼt = max(ℯ * abs(agepointsₚ[k] - agepoints[k]), dt)
                    end
                    if Tpointsₚ[k] != Tpoints[k]
                        σⱼT = max(ℯ * abs(Tpointsₚ[k] - Tpoints[k]), one(T))
                    end
                end

                # Update the currently accepted proposal
                ll = llₚ
                npoints = npointsₚ
                copyto!(agepoints, agepointsₚ)
                copyto!(Tpoints, Tpointsₚ)
                copyto!(unconf.agepoints, unconf.agepointsₚ)
                copyto!(unconf.Tpoints, unconf.Tpointsₚ)
                copyto!(boundary.Tpoints, boundary.Tpointsₚ)
                copyto!(calcHeAges, calcHeAgesₚ)

                # Not critical to the function of the MCMC loop, but critical for recording stationary distribution!
                copyto!(Tsteps, Tstepsₚ)
                acceptancedist[n] = true
            end

            # Record results for analysis and troubleshooting
            lldist[n] = normpdf_ll(HeAge, σ, calcHeAges) # Recalculated to constant baseline
            ndist[n] = npoints # distribution of # of points
            σⱼtdist[n] = σⱼt
            σⱼTdist[n] = σⱼT
            HeAgedist[:,n] .= calcHeAges # distribution of He ages

            # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
            collectto!(view(tpointdist, :, n), view(agepoints, 1:npoints), boundary.agepoints, unconf.agepoints)
            collectto!(view(Tpointdist, :, n), view(Tpoints, 1:npoints), boundary.Tpoints, unconf.Tpoints)

            # Update progress meter every `progress_interval` steps
            (mod(n, progress_interval) == 0) && update!(progress, n)
        end
        return (tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist, σⱼtdist, σⱼTdist)
    end
    export MCMC_vartcryst


    function MCMC_vartcryst(data::NamedTuple, model::NamedTuple, npoints::Int, agepoints::DenseVector{T}, Tpoints::DenseVector{T}, boundary::Boundary{T}, unconf::Unconformity{T}, detail::DetailInterval{T}) where T <: AbstractFloat
        # Sanitize inputs
        @assert firstindex(agepoints) === 1
        @assert firstindex(Tpoints) === 1
        halfwidth = T.(data.halfwidth)::DenseVector{T}
        halfwidth = T.(data.halfwidth)::DenseVector{T}
        crystAge = T.(data.crystAge)::DenseVector{T}
        HeAge = T.(data.HeAge)::DenseVector{T}
        HeAge_sigma = T.(data.HeAge_sigma)::DenseVector{T}
        U = T.(data.U)::DenseVector{T}
        Th = T.(data.Th)::DenseVector{T}
        nsteps = model.nsteps::Int
        maxpoints = model.maxpoints::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        totalpoints = maxpoints + boundary.npoints + unconf.npoints::Int
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(model.dTmax)::T
        agesteps = T.(model.agesteps)::DenseVector{T}
        tsteps = T.(model.tsteps)::DenseVector{T}
        tinit = T(model.tinit)::T
        Tinit = T(model.Tinit)::T
        Tnow = T(model.Tnow)::T
        dt = T(model.dt)::T
        dr = T(model.dr)::T
        σmodel = T(model.σmodel)::T
        σannealing = T(model.σannealing)::T
        λannealing = T(model.λannealing)::T

        # Calculate number of boundary and unconformity points and allocate buffer for interpolating
        agepointbuffer = similar(agepoints, totalpoints)::DenseVector{T}
        Tpointbuffer = similar(agepoints, totalpoints)::DenseVector{T}
        knot_index = similar(agesteps, Int)::DenseVector{Int}

        # Calculate model ages for initial proposal
        ages = collectto!(agepointbuffer, view(agepoints, 1:npoints), boundary.agepoints, unconf.agepoints)::StridedVector{T}
        temperatures = collectto!(Tpointbuffer, view(Tpoints, 1:npoints), boundary.Tpoints, unconf.Tpoints)::StridedVector{T}
        Tsteps = linterp1s(ages, temperatures, agesteps)::DenseVector{T}
        calcHeAges = similar(HeAge)::DenseVector{T}
        pr, Teq = anneal(dt, tsteps, Tsteps, ZRDAAM()) # Damage annealing history
        pr::DenseMatrix{T}
        Teq::DenseVector{T}

        # Prepare a Mineral object for each analysis
        diffusionmodel = (DzEa=T(model.DzEa), DzD0=T(model.DzD0), DN17Ea=T(model.DN17Ea), DN17D0=T(model.DN17D0))
        zircons = Array{Zircon{T}}(undef, length(halfwidth))::Vector{Zircon{T}}
        for i=1:length(zircons)
            # Iterate through each grain, calculate the modeled age for each
            first_index = 1 + floor(Int64,(tinit - crystAge[i])/dt)
            zircons[i] = Zircon(halfwidth[i], dr, U[i], Th[i], dt, agesteps[first_index:end])
            calcHeAges[i] = HeAgeSpherical(zircons[i], @views(Tsteps[first_index:end]), @views(pr[first_index:end,first_index:end]), diffusionmodel)::T
        end

        # Simulated annealing of uncertainty
        σₐ = simannealsigma.(1, HeAge_sigma, σmodel, σannealing, λannealing)::DenseVector{T}
        σ = sqrt.(HeAge_sigma.^2 .+ σmodel^2)

        # Log-likelihood for initial proposal
        ll = normpdf_ll(HeAge, σₐ, calcHeAges)
        if simplified
            ll -= log(npoints)
        end

        # Variables to hold proposals
        llₚ = ll
        npointsₚ = npoints
        agepointsₚ = copy(agepoints)::DenseVector{T}
        Tpointsₚ = copy(Tpoints)::DenseVector{T}
        calcHeAgesₚ = copy(calcHeAges)::DenseVector{T}
        Tstepsₚ = copy(Tsteps)::DenseVector{T}

        # distributions to populate
        tpointdist = zeros(T, totalpoints, nsteps)
        Tpointdist = zeros(T, totalpoints, nsteps)
        HeAgedist = zeros(T, length(HeAge), nsteps)
        σⱼtdist = zeros(T, nsteps)
        σⱼTdist = zeros(T, nsteps)
        lldist = zeros(T, nsteps)
        ndist = zeros(Int, nsteps)
        acceptancedist = zeros(Bool, nsteps)

        # Standard deviations of Gaussian proposal ("jumping") distributions
        # for temperature and time
        σⱼt = tinit/60
        σⱼT = Tinit/60
        k = 1 # Index of chosen t-T point

        # Proposal probabilities (must sum to 1)
        move = 0.64
        birth = 0.15 # Must equal death
        death = 0.15 # Must equal birth
        movebounds = 0.06

        # Number of times to attempt to satisfy reheating rate: # Should be large
        # enough that proposal probabilities are unchanged, but low enough to prevent
        # an infinite loop
        nattempts = 100_000

        progress = Progress(nsteps, dt=1, desc="Running MCMC ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps

            # Copy proposal from last accepted solution
            npointsₚ = npoints
            copyto!(agepointsₚ, agepoints)
            copyto!(Tpointsₚ, Tpoints)
            copyto!(unconf.agepointsₚ, unconf.agepoints)
            copyto!(unconf.Tpointsₚ, unconf.Tpoints)
            copyto!(boundary.Tpointsₚ, boundary.Tpoints)
            enoughpoints = min(pointsininterval(agepoints, npoints, detail.agemin, detail.agemax), detail.minpoints)::Int

            # Adjust the proposal
            r = rand()
            if r < move
                # Move one t-T point
                k = ceil(Int, rand() * npoints)
                for attempt ∈ 1:nattempts

                # Move the age of one model point
                agepointsₚ[k] += randn() * σⱼt
                # if agepointsₚ[k] < (0 + dt)
                #     # Reflecting boundary condition at (0 + dt)
                #     agepointsₚ[k] = (0 + dt) - (agepointsₚ[k] - (0 + dt))
                # elseif agepointsₚ[k] > (tinit - dt)
                #     # Reflecting boundary condition at (tinit - dt)
                #     agepointsₚ[k] = (tinit - dt) - (agepointsₚ[k] - (tinit - dt))
                # end

                # Move the Temperature of one model point
                Tpointsₚ[k] += randn() * σⱼT
                # if Tpointsₚ[k] < Tnow
                #     # Reflecting boundary conditions at Tnow (0)
                #     Tpointsₚ[k] = Tnow - (Tpointsₚ[k] - Tnow)
                # elseif Tpointsₚ[k] > Tinit
                #     # Reflecting boundary conditions at Tinit
                #     Tpointsₚ[k] = Tinit - (Tpointsₚ[k] - Tinit)
                # end

                # Circular boundary conditions
                agepointsₚ[k] = mod(agepointsₚ[k]-dt, tinit-2dt) + dt
                Tpointsₚ[k] = mod(Tpointsₚ[k]-Tnow, Tinit-Tnow) + Tnow
                # agepointsₚ[k] = min(max(agepointsₚ[k], dt), tinit)
                # Tpointsₚ[k] = min(max(Tpointsₚ[k], Tnow), Tinit)

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, unconf.agepointsₚ)::StridedVector{T}
                temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, unconf.Tpointsₚ)::StridedVector{T}
                linterp1s!(Tstepsₚ, knot_index, ages, temperatures, agesteps)

                # Retry unless we have satisfied the maximum reheating rate
                if isdistinct(agepointsₚ, npointsₚ, k, 2dt) && maxdiff(Tstepsₚ) < dTmax
                    if pointsininterval(agepointsₚ, npointsₚ, detail.agemin, detail.agemax) >= enoughpoints
                        break
                    end
                end
                if attempt == nattempts
                    @warn """`move` proposals failed to satisfy reheating rate limit
                    maxdiff: $(maxdiff(Tstepsₚ))
                    ages: $(agepointsₚ[1:npointsₚ])
                    temperatures: $(Tpointsₚ[1:npointsₚ])
                    σⱼt: $(σⱼt)
                    σⱼT: $(σⱼT)"""
                end
                # Copy last accepted solution to re-modify if we don't break
                copyto!(agepointsₚ, agepoints)
                copyto!(Tpointsₚ, Tpoints)
                end
            elseif (r < move+birth) && (npointsₚ < maxpoints)
                # Birth: add a new model point
                npointsₚ = npoints + 1
                for attempt ∈ 1:nattempts
                agepointsₚ[npointsₚ] = 0 + rand()*(tinit-0)
                Tpointsₚ[npointsₚ] = Tnow + rand()*(Tinit-Tnow)

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, unconf.agepointsₚ)::StridedVector{T}
                temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, unconf.Tpointsₚ)::StridedVector{T}
                linterp1s!(Tstepsₚ, knot_index, ages, temperatures, agesteps)

                # Retry unless we have satisfied the maximum reheating rate
                if isdistinct(agepointsₚ, npointsₚ, npointsₚ, 2dt) && maxdiff(Tstepsₚ) < dTmax
                    break
                end
                if attempt == nattempts
                    @warn """new point proposals failed to satisfy reheating rate limit
                    maxdiff: $(maxdiff(Tstepsₚ))
                    ages: $(agepointsₚ[1:npointsₚ])
                    temperatures: $(Tpointsₚ[1:npointsₚ])
                    σⱼt: $(σⱼt)
                    σⱼT: $(σⱼT)"""
                end
                end
            elseif (r < move+birth+death) && (r >= move+birth) && (npoints > max(minpoints, detail.minpoints))
                # Death: remove a model point
                npointsₚ = npoints - 1 # Delete last point in array from proposal
                for attempt ∈ 1:nattempts
                k = ceil(Int, rand()*npoints) # Choose point to delete
                agepointsₚ[k] = agepointsₚ[npoints]
                Tpointsₚ[k] = Tpointsₚ[npoints]

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, unconf.agepointsₚ)::StridedVector{T}
                temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, unconf.Tpointsₚ)::StridedVector{T}
                linterp1s!(Tstepsₚ, knot_index, ages, temperatures, agesteps)

                # Retry unless we have satisfied the maximum reheating rate
                if maxdiff(Tstepsₚ) < dTmax
                    if pointsininterval(agepointsₚ, npointsₚ, detail.agemin, detail.agemax) >= enoughpoints
                        break
                    end
                end
                if attempt == nattempts
                    @warn """point removal proposals failed to satisfy reheating rate limit
                    maxdiff: $(maxdiff(Tstepsₚ))
                    ages: $(agepointsₚ[1:npointsₚ])
                    temperatures: $(Tpointsₚ[1:npointsₚ])
                    σⱼt: $(σⱼt)
                    σⱼT: $(σⱼT)"""
                end
                # Copy last accepted solution to re-modify if we don't break
                copyto!(agepointsₚ, agepoints)
                copyto!(Tpointsₚ, Tpoints)
                end
            else
                # Move boundary conditions
                for attempt ∈ 1:nattempts
                # Move the temperatures of the starting or ending boundaries
                # @. boundary.Tpointsₚ = boundary.T₀ + rand()*boundary.ΔT
                k = rand(1:boundary.npoints)
                boundary.Tpointsₚ[k] = boundary.T₀[k] + rand()*boundary.ΔT[k]

                # If there's an imposed unconformity, adjust within parameters
                if unconf.npoints > 0
                    @. unconf.agepointsₚ = unconf.Age₀ + rand()*unconf.ΔAge
                    @. unconf.Tpointsₚ = unconf.T₀ + rand()*unconf.ΔT
                end

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, unconf.agepointsₚ)::StridedVector{T}
                temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, unconf.Tpointsₚ)::StridedVector{T}
                linterp1s!(Tstepsₚ, knot_index, ages, temperatures, agesteps)

                # Retry unless we have satisfied the maximum reheating rate
                (maxdiff(Tstepsₚ) < dTmax) && break

                # Copy last accepted solution to re-modify if we don't break
                if attempt == nattempts
                    @warn """`movebounds` proposals failed to satisfy reheating rate limit
                    maxdiff: $(maxdiff(Tstepsₚ))
                    ages: $(agepointsₚ[1:npointsₚ])
                    temperatures: $(Tpointsₚ[1:npointsₚ])
                    σⱼt: $(σⱼt)
                    σⱼT: $(σⱼT)"""
                end
                copyto!(unconf.agepointsₚ, unconf.agepoints)
                copyto!(unconf.Tpointsₚ, unconf.Tpoints)
                copyto!(boundary.Tpointsₚ, boundary.Tpoints)
                end
            end

            if any(isnan, view(agepointsₚ, 1:npointsₚ)) ||  any(isnan, view(Tpointsₚ, 1:npointsₚ))
                @warn """`NaN`s detected!
                ages: $(agepointsₚ[1:npointsₚ])
                temperatures: $(Tpointsₚ[1:npointsₚ])
                σⱼt: $(σⱼt)
                σⱼT: $(σⱼT)
                r: $r
                """
            end

             # Calculate model ages for each grain
            anneal!(pr, Teq, dt, tsteps, Tstepsₚ, ZRDAAM())
            for i=1:length(zircons)
                first_index = 1 + floor(Int64,(tinit - crystAge[i])/dt)
                calcHeAgesₚ[i] = HeAgeSpherical(zircons[i], @views(Tstepsₚ[first_index:end]), @views(pr[first_index:end,first_index:end]), diffusionmodel)::T
            end

            # Calculate log likelihood of proposal
            σₐ .= simannealsigma.(n, HeAge_sigma, σmodel, σannealing, λannealing)
            llₚ = normpdf_ll(HeAge, σₐ, calcHeAgesₚ)
            llₗ = normpdf_ll(HeAge, σₐ, calcHeAges) # Recalulate last one too with new σₐ
            if simplified # slightly penalize more complex t-T paths
                llₚ -= log(npointsₚ)
                llₗ -= log(npoints)
            end

            # Accept or reject proposal based on likelihood
            # To avoid numerical problems with diffusion code, also reject proposal
            # if maximum proposed heating rate is greater than 25C per timestep.
            # (Fast cooling should not be a problem, however)
            if log(rand()) < (llₚ - llₗ)

                # Update jumping distribution based on size of current accepted move
                if dynamicjumping && r < move
                    if agepointsₚ[k] != agepoints[k]
                        σⱼt = max(ℯ * abs(agepointsₚ[k] - agepoints[k]), dt)
                    end
                    if Tpointsₚ[k] != Tpoints[k]
                        σⱼT = max(ℯ * abs(Tpointsₚ[k] - Tpoints[k]), one(T))
                    end
                end

                # Update the currently accepted proposal
                ll = llₚ
                npoints = npointsₚ
                copyto!(agepoints, agepointsₚ)
                copyto!(Tpoints, Tpointsₚ)
                copyto!(unconf.agepoints, unconf.agepointsₚ)
                copyto!(unconf.Tpoints, unconf.Tpointsₚ)
                copyto!(boundary.Tpoints, boundary.Tpointsₚ)
                copyto!(calcHeAges, calcHeAgesₚ)

                # Not critical to the function of the MCMC loop, but critical for recording stationary distribution!
                copyto!(Tsteps, Tstepsₚ)
                acceptancedist[n] = true
            end

            # Record results for analysis and troubleshooting
            lldist[n] = normpdf_ll(HeAge, σ, calcHeAges) # Recalculated to constant baseline
            ndist[n] = npoints # distribution of # of points
            σⱼtdist[n] = σⱼt
            σⱼTdist[n] = σⱼT
            HeAgedist[:,n] .= calcHeAges # distribution of He ages

            # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
            collectto!(view(tpointdist, :, n), view(agepoints, 1:npoints), boundary.agepoints, unconf.agepoints)
            collectto!(view(Tpointdist, :, n), view(Tpoints, 1:npoints), boundary.Tpoints, unconf.Tpoints)

            # Update progress meter every `progress_interval` steps
            (mod(n, progress_interval) == 0) && update!(progress, n)
        end
        return (tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist, σⱼtdist, σⱼTdist)
    end

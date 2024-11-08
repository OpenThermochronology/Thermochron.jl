
    """
    ```julia
    MCMC(data::NamedTuple, model::NamedTuple, npoints::Int, agepoints::Vector, Tpoints::Vector, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])
    ```
    Markov chain Monte Carlo time-Temperature inversion of the data specified in `data` and model parameters specified by `model`.

    Returns a tuple of distributions `(tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist)`

    ## Examples
    ```julia
    tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist = MCMC(data, model, npoints, agepoints, Tpoints, constraint, boundary)
    ```
    """
    function MCMC(data::NamedTuple, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where T <: AbstractFloat
        # Process inputs
        halfwidth = T.(data.halfwidth)::DenseVector{T}
        crystAge = T.(data.crystAge)::DenseVector{T}
        HeAge = T.(data.HeAge)::DenseVector{T}
        HeAge_sigma = T.(data.HeAge_sigma)::DenseVector{T}
        U = T.(data.U)::DenseVector{T}
        Th = T.(data.Th)::DenseVector{T}
        Sm = (haskey(data, :Sm) ? T.(data.Th) : zeros(T, size(U)))::DenseVector{T}
        nsteps = (haskey(model, :nsteps) ? model.nsteps : 10^6)::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        maxpoints = (haskey(model, :maxpoints) ? model.maxpoints : 50)::Int
        npoints = (haskey(model, :npoints) ? model.npoints : minpoints)::Int
        totalpoints = maxpoints + boundary.npoints + constraint.npoints::Int
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : 5)::T
        agesteps = T.(model.agesteps)::DenseVector{T}
        tsteps = T.(model.tsteps)::DenseVector{T}
        tnow, tinit = extrema(boundary.agepoints)
        Tnow, Tinit = extrema(boundary.T₀)
        Tr = T(haskey(model, :Tr) ? model.Tr : (Tinit+Tnow)/2)::T
        dt = T(model.dt)::T
        dr = T(model.dr)::T
        σmodel = T(model.σmodel)::T
        σannealing = T(model.σannealing)::T
        λannealing = T(model.λannealing)::T

        # Arrays to hold all t and T points (up to npoints=maxpoints)
        agepoints = zeros(T, maxpoints) 
        Tpoints = zeros(T, maxpoints)
    
        # Fill some intermediate points to give the MCMC something to work with
        agepoints[1:npoints] .= range(tnow, tinit, length=npoints)
        Tpoints[1:npoints] .= Tr # Degrees C

        # Calculate number of boundary and unconformity points and allocate buffer for interpolating
        agepointbuffer = similar(agepoints, totalpoints)::DenseVector{T}
        Tpointbuffer = similar(agepoints, totalpoints)::DenseVector{T}
        knot_index = similar(agesteps, Int)::DenseVector{Int}

        # Calculate model ages for initial proposal
        ages = collectto!(agepointbuffer, view(agepoints, 1:npoints), boundary.agepoints, constraint.agepoints)::StridedVector{T}
        temperatures = collectto!(Tpointbuffer, view(Tpoints, 1:npoints), boundary.Tpoints, constraint.Tpoints)::StridedVector{T}
        Tsteps = linterp1s(ages, temperatures, agesteps)::DenseVector{T}
        calcHeAges = similar(HeAge)::DenseVector{T}

        # Damage models for each mineral
        zdm = if haskey(model, :DzEa) && haskey(model, :DzD0) && haskey(model, :DN17Ea) && haskey(model, :DN17D0)
            ZRDAAM(
                DzEa=T(model.DzEa), 
                DzD0=T(model.DzD0), 
                DN17Ea=T(model.DN17Ea), 
                DN17D0=T(model.DN17D0)
            )
        else
            ZRDAAM()
        end
        zpr, zTeq = anneal(dt, tsteps, Tsteps, zdm) # Zircon amage annealing history
        zpr::DenseMatrix{T}
        zTeq::DenseVector{T}

        adm = RDAAM()
        apr, aTeq = anneal(dt, tsteps, Tsteps, adm) # Apatite damage annealing history
        apr::DenseMatrix{T}
        aTeq::DenseVector{T}

        # See what minerals we have
        tzr = containsi.(data.mineral, "zircon")
        any(tzr) && @info "Inverting for He ages of $(count(tzr)) zircons"
        tap = containsi.(data.mineral, "apatite")
        any(tap) && @info "Inverting for He ages of $(count(tap)) apatites"

        # Prepare a Mineral object for each analysis
        zircons = Array{Zircon{T}}(undef, count(tzr))::Vector{Zircon{T}}
        zi = 1
        for i in findall(tzr)
            # Iterate through each grain, calculate the modeled age for each
            first_index = 1 + floor(Int64,(tinit - crystAge[i])/dt)
            zircons[zi] = Zircon(halfwidth[i], dr, U[i], Th[i], Sm[i], dt, agesteps[first_index:end])
            calcHeAges[i] = HeAgeSpherical(zircons[zi], @views(Tsteps[first_index:end]), @views(zpr[first_index:end,first_index:end]), zdm)::T
            zi += 1
        end
        apatites = Array{Apatite{T}}(undef, count(tap))::Vector{Apatite{T}}
        ai = 1
        for i in findall(tap)
            # Iterate through each grain, calculate the modeled age for each
            first_index = 1 + floor(Int64,(tinit - crystAge[i])/dt)
            apatites[ai] = Apatite(halfwidth[i], dr, U[i], Th[i], Sm[i], dt, agesteps[first_index:end])
            calcHeAges[i] = HeAgeSpherical(apatites[ai], @views(Tsteps[first_index:end]), @views(apr[first_index:end,first_index:end]), adm)::T
            ai += 1
        end

        # Simulated annealing of uncertainty
        σₐ = simannealsigma.(1, HeAge_sigma, σmodel, σannealing, λannealing)::DenseVector{T}
        σ = sqrt.(HeAge_sigma.^2 .+ σmodel^2)

        # Log-likelihood for initial proposal
        llna = llnaₚ = diff_ll(Tsteps, dTmax, dTmax_sigma) + (simplified ? -log(npoints) : zero(T))
        ll = llₚ =  normpdf_ll(HeAge, σₐ, calcHeAges) + llna

        # Variables to hold proposals
        npointsₚ = npoints
        agepointsₚ = copy(agepoints)::DenseVector{T}
        Tpointsₚ = copy(Tpoints)::DenseVector{T}
        calcHeAgesₚ = copy(calcHeAges)::DenseVector{T}

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
        σⱼt = fill((tinit-tnow)/60, maxpoints)
        σⱼT = fill((Tinit-Tnow)/60, maxpoints)
        k = 1 # Index of chosen t-T point

        # Proposal probabilities (must sum to 1)
        p_move = 0.64
        p_birth = 0.15 # Must equal p_death
        p_death = 0.15 # Must equal p_birth
        p_bounds = 0.06

        progress = Progress(nsteps, dt=1, desc="Running MCMC ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps
            if detail.minpoints > 0
                enoughpoints = min(pointsininterval(agepoints, npoints, detail.agemin, detail.agemax, dt), detail.minpoints)::Int
            end
            @label restart

            # Copy proposal from last accepted solution
            npointsₚ = npoints
            copyto!(agepointsₚ, agepoints)
            copyto!(Tpointsₚ, Tpoints)
            copyto!(constraint.agepointsₚ, constraint.agepoints)
            copyto!(constraint.Tpointsₚ, constraint.Tpoints)
            copyto!(boundary.Tpointsₚ, boundary.Tpoints)

            # Randomly choose an option and point (if applicable) to adjust
            r = rand()
            k = ceil(Int, rand()*npoints) 

            # Adjust the proposal
            if r < p_move
                # Move one t-T point
                movepoint!(agepointsₚ, Tpointsₚ, k, σⱼt[k], σⱼT[k], boundary)

            elseif (r < p_move+p_birth) && (npoints < maxpoints)
                # Birth: add a new model point
                k = npointsₚ = npoints + 1
                agepointsₚ[k] = tnow + rand()*(tinit-tnow)
                Tpointsₚ[k] = Tnow + rand()*(Tinit-Tnow)
                σⱼt[k] = (tinit-Tnow)/60
                σⱼT[k] = (Tinit-Tnow)/60

            elseif (r < p_move+p_birth+p_death) && (r >= p_move+p_birth) && (npoints > max(minpoints, detail.minpoints))
                # Death: remove a model point
                npointsₚ = npoints - 1
                agepointsₚ[k] = agepointsₚ[npoints]
                Tpointsₚ[k] = Tpointsₚ[npoints]

            elseif (r < p_move+p_birth+p_death+p_bounds)
                # Move the temperatures of the starting and ending boundaries
                movebounds!(boundary)
                # If there's an imposed unconformity or other t-T constraint, adjust within parameters
                movebounds!(constraint)

            end

            # Recalculate interpolated proposed t-T path
            ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, constraint.agepointsₚ)::StridedVector{T}
            temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, constraint.Tpointsₚ)::StridedVector{T}
            linterp1s!(Tsteps, knot_index, ages, temperatures, agesteps)
            
            # Old version of imposing max reheating rate, for reference:
            # (maxdiff(Tsteps) > dTmax) && @goto restart

            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (pointsininterval(agepointsₚ, npointsₚ, detail.agemin, detail.agemax, dt) < enoughpoints) && @goto restart
            end
            # Calculate model ages for each grain
            if any(tzr)
                anneal!(zpr, zTeq, dt, tsteps, Tsteps, zdm)
                zi = 1
                for i ∈ findall(tzr)
                    first_index = 1 + floor(Int64,(tinit - crystAge[i])/dt)
                    calcHeAgesₚ[i] = HeAgeSpherical(zircons[zi], @views(Tsteps[first_index:end]), @views(zpr[first_index:end,first_index:end]), zdm)::T
                    zi += 1
                end
            end
            if any(tap)
                anneal!(apr, aTeq, dt, tsteps, Tsteps, adm)
                ai = 1
                for i ∈ findall(tap)
                    first_index = 1 + floor(Int64,(tinit - crystAge[i])/dt)
                    calcHeAgesₚ[i] = HeAgeSpherical(apatites[ai], @views(Tsteps[first_index:end]), @views(apr[first_index:end,first_index:end]), adm)::T
                    ai += 1
                end
            end

            # Calculate log likelihood of proposal
            llnaₚ = diff_ll(Tsteps, dTmax, dTmax_sigma)
            simplified && (llnaₚ += -log(npointsₚ))
            σₐ .= simannealsigma.(n, HeAge_sigma, σmodel, σannealing, λannealing)
            llₚ = normpdf_ll(HeAge, σₐ, calcHeAgesₚ) + llnaₚ
            llₗ = normpdf_ll(HeAge, σₐ, calcHeAges) + llna # Recalulate last one too with new σₐ

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - llₗ)

                # Update jumping distribution based on size of current accepted p_move
                if dynamicjumping && r < p_move
                    if agepointsₚ[k] != agepoints[k]
                        σⱼt[k] = ℯ * abs(agepointsₚ[k] - agepoints[k])
                    end
                    if Tpointsₚ[k] != Tpoints[k]
                        σⱼT[k] = ℯ * abs(Tpointsₚ[k] - Tpoints[k])
                    end
                end

                # Update the currently accepted proposal
                ll = llₚ
                llna = llnaₚ
                npoints = npointsₚ
                copyto!(agepoints, 1, agepointsₚ, 1, npoints)
                copyto!(Tpoints, 1, Tpointsₚ, 1, npoints)
                copyto!(constraint.agepoints, constraint.agepointsₚ)
                copyto!(constraint.Tpoints, constraint.Tpointsₚ)
                copyto!(boundary.Tpoints, boundary.Tpointsₚ)
                copyto!(calcHeAges, calcHeAgesₚ)

                # Not critical to the function of the MCMC loop, but critical for recording stationary distribution!
                acceptancedist[n] = true
            end

            # Record results for analysis and troubleshooting
            lldist[n] = llna + normpdf_ll(HeAge, σ, calcHeAges) # Recalculated to constant baseline
            ndist[n] = npoints # distribution of # of points
            σⱼtdist[n] = σⱼt[k]
            σⱼTdist[n] = σⱼT[k]
            HeAgedist[:,n] .= calcHeAges # distribution of He ages

            # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
            collectto!(view(tpointdist, :, n), view(agepoints, 1:npoints), boundary.agepoints, constraint.agepoints)
            collectto!(view(Tpointdist, :, n), view(Tpoints, 1:npoints), boundary.Tpoints, constraint.Tpoints)

            # Update progress meter every `progress_interval` steps
            (mod(n, progress_interval) == 0) && update!(progress, n)
        end
        return (tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist, σⱼtdist, σⱼTdist)
    end
    export MCMC

    function MCMC_varkinetics(data::NamedTuple, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where T <: AbstractFloat
        # Process inputs
        halfwidth = T.(data.halfwidth)::DenseVector{T}
        crystAge = T.(data.crystAge)::DenseVector{T}
        HeAge = T.(data.HeAge)::DenseVector{T}
        HeAge_sigma = T.(data.HeAge_sigma)::DenseVector{T}
        U = T.(data.U)::DenseVector{T}
        Th = T.(data.Th)::DenseVector{T}
        Sm = (haskey(data, :Sm) ? T.(data.Th) : zeros(T, size(U)))::DenseVector{T}
        nsteps = (haskey(model, :nsteps) ? model.nsteps : 10^6)::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        maxpoints = (haskey(model, :maxpoints) ? model.maxpoints : 50)::Int
        npoints = (haskey(model, :npoints) ? model.npoints : minpoints)::Int
        totalpoints = maxpoints + boundary.npoints + constraint.npoints::Int
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : 5)::T
        agesteps = T.(model.agesteps)::DenseVector{T}
        tsteps = T.(model.tsteps)::DenseVector{T}
        tnow, tinit = extrema(boundary.agepoints)
        Tnow, Tinit = extrema(boundary.T₀)
        Tr = T(haskey(model, :Tr) ? model.Tr : (Tinit+Tnow)/2)::T
        dt = T(model.dt)::T
        dr = T(model.dr)::T
        σmodel = T(model.σmodel)::T
        σannealing = T(model.σannealing)::T
        λannealing = T(model.λannealing)::T

        # Arrays to hold all t and T points (up to npoints=maxpoints)
        agepoints = zeros(T, maxpoints) 
        Tpoints = zeros(T, maxpoints)
    
        # Fill some intermediate points to give the MCMC something to work with
        agepoints[1:npoints] .= range(tnow, tinit, length=npoints)
        Tpoints[1:npoints] .= Tr # Degrees C

        # Calculate number of boundary and unconformity points and allocate buffer for interpolating
        agepointbuffer = similar(agepoints, totalpoints)::DenseVector{T}
        Tpointbuffer = similar(agepoints, totalpoints)::DenseVector{T}
        knot_index = similar(agesteps, Int)::DenseVector{Int}

        # Calculate model ages for initial proposal
        ages = collectto!(agepointbuffer, view(agepoints, 1:npoints), boundary.agepoints, constraint.agepoints)::StridedVector{T}
        temperatures = collectto!(Tpointbuffer, view(Tpoints, 1:npoints), boundary.Tpoints, constraint.Tpoints)::StridedVector{T}
        Tsteps = linterp1s(ages, temperatures, agesteps)::DenseVector{T}
        calcHeAges = similar(HeAge)::DenseVector{T}

        # Damage models for each mineral
        zdm₀ = zdm = zdmₚ = if haskey(model, :DzEa) && haskey(model, :DzD0) && haskey(model, :DN17Ea) && haskey(model, :DN17D0)
            ZRDAAM(
                DzEa=T(model.DzEa), 
                DzD0=T(model.DzD0), 
                DN17Ea=T(model.DN17Ea), 
                DN17D0=T(model.DN17D0)
            )
        else
            ZRDAAM()
        end
        zpr, zTeq = anneal(dt, tsteps, Tsteps, zdm) # Zircon amage annealing history
        zpr::DenseMatrix{T}
        zTeq::DenseVector{T}

        adm₀ = adm = admₚ = RDAAM()
        apr, aTeq = anneal(dt, tsteps, Tsteps, adm) # Apatite damage annealing history
        apr::DenseMatrix{T}
        aTeq::DenseVector{T}

        # See what minerals we have
        tzr = containsi.(data.mineral, "zircon")
        any(tzr) && @info "Inverting for He ages of $(count(tzr)) zircons"
        tap = containsi.(data.mineral, "apatite")
        any(tap) && @info "Inverting for He ages of $(count(tap)) apatites"

        # Prepare a Mineral object for each analysis
        zircons = Array{Zircon{T}}(undef, count(tzr))::Vector{Zircon{T}}
        zi = 1
        for i in findall(tzr)
            # Iterate through each grain, calculate the modeled age for each
            first_index = 1 + floor(Int64,(tinit - crystAge[i])/dt)
            zircons[zi] = Zircon(halfwidth[i], dr, U[i], Th[i], Sm[i], dt, agesteps[first_index:end])
            calcHeAges[i] = HeAgeSpherical(zircons[zi], @views(Tsteps[first_index:end]), @views(zpr[first_index:end,first_index:end]), zdm)::T
            zi += 1
        end
        apatites = Array{Apatite{T}}(undef, count(tap))::Vector{Apatite{T}}
        ai = 1
        for i in findall(tap)
            # Iterate through each grain, calculate the modeled age for each
            first_index = 1 + floor(Int64,(tinit - crystAge[i])/dt)
            apatites[ai] = Apatite(halfwidth[i], dr, U[i], Th[i], Sm[i], dt, agesteps[first_index:end])
            calcHeAges[i] = HeAgeSpherical(apatites[ai], @views(Tsteps[first_index:end]), @views(apr[first_index:end,first_index:end]), adm)::T
            ai += 1
        end

        # Simulated annealing of uncertainty
        σₐ = simannealsigma.(1, HeAge_sigma, σmodel, σannealing, λannealing)::DenseVector{T}
        σ = sqrt.(HeAge_sigma.^2 .+ σmodel^2)

        # Log-likelihood for initial proposal
        llna = llnaₚ = diff_ll(Tsteps, dTmax, dTmax_sigma) + loglikelihood(admₚ, adm₀) + loglikelihood(zdmₚ, zdm₀) + (simplified ? -log(npoints) : zero(T))
        ll = llₚ =  normpdf_ll(HeAge, σₐ, calcHeAges) + llna

        # Variables to hold proposals
        npointsₚ = npoints
        agepointsₚ = copy(agepoints)::DenseVector{T}
        Tpointsₚ = copy(Tpoints)::DenseVector{T}
        calcHeAgesₚ = copy(calcHeAges)::DenseVector{T}

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
        σⱼt = fill((tinit-tnow)/60, maxpoints)
        σⱼT = fill((Tinit-Tnow)/60, maxpoints)
        k = 1 # Index of chosen t-T point

        # Proposal probabilities (must sum to 1)
        p_move = 0.6
        p_birth = 0.15 # Must equal p_death
        p_death = 0.15 # Must equal p_birth
        p_bounds = 0.06
        p_kinetics = 0.04

        progress = Progress(nsteps, dt=1, desc="Running MCMC ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps
            if detail.minpoints > 0
                enoughpoints = min(pointsininterval(agepoints, npoints, detail.agemin, detail.agemax, dt), detail.minpoints)::Int
            end
            @label restart

            # Copy proposal from last accepted solution
            admₚ = adm
            zdmₚ = zdm
            npointsₚ = npoints
            copyto!(agepointsₚ, agepoints)
            copyto!(Tpointsₚ, Tpoints)
            copyto!(constraint.agepointsₚ, constraint.agepoints)
            copyto!(constraint.Tpointsₚ, constraint.Tpoints)
            copyto!(boundary.Tpointsₚ, boundary.Tpoints)

            # Randomly choose an option and point (if applicable) to adjust
            r = rand()
            k = ceil(Int, rand()*npoints) 

            # Adjust the proposal
            if r < p_move
                # Move one t-T point
                movepoint!(agepointsₚ, Tpointsₚ, k, σⱼt[k], σⱼT[k], boundary)

            elseif (r < p_move+p_birth) && (npoints < maxpoints)
                # Birth: add a new model point
                k = npointsₚ = npoints + 1
                agepointsₚ[k] = tnow + rand()*(tinit-tnow)
                Tpointsₚ[k] = Tnow + rand()*(Tinit-Tnow)
                σⱼt[k] = (tinit-Tnow)/60
                σⱼT[k] = (Tinit-Tnow)/60

            elseif (r < p_move+p_birth+p_death) && (r >= p_move+p_birth) && (npoints > max(minpoints, detail.minpoints))
                # Death: remove a model point
                npointsₚ = npoints - 1
                agepointsₚ[k] = agepointsₚ[npoints]
                Tpointsₚ[k] = Tpointsₚ[npoints]

            elseif (r < p_move+p_birth+p_death+p_bounds)
                # Move the temperatures of the starting and ending boundaries
                movebounds!(boundary)
                # If there's an imposed unconformity or other t-T constraint, adjust within parameters
                movebounds!(constraint)

            elseif (r < p_move+p_birth+p_death+p_bounds+p_kinetics)
                # Adjust kinetic parameters, one at a time
                any(tzr) && (zdmₚ = movekinetics(zdm))
                any(tap) && (admₚ = movekinetics(adm))

            end

            # Recalculate interpolated proposed t-T path
            ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, constraint.agepointsₚ)::StridedVector{T}
            temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, constraint.Tpointsₚ)::StridedVector{T}
            linterp1s!(Tsteps, knot_index, ages, temperatures, agesteps)
    
            # Old version of imposing max reheating rate, for reference:
            # (maxdiff(Tsteps) > dTmax) && @goto restart 

            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (pointsininterval(agepointsₚ, npointsₚ, detail.agemin, detail.agemax, dt) < enoughpoints) && @goto restart
            end
               
            # Calculate model ages for each grain
            if any(tzr)
                anneal!(zpr, zTeq, dt, tsteps, Tsteps, zdmₚ)
                zi = 1
                for i ∈ findall(tzr)
                    first_index = 1 + floor(Int64,(tinit - crystAge[i])/dt)
                    calcHeAgesₚ[i] = HeAgeSpherical(zircons[zi], @views(Tsteps[first_index:end]), @views(zpr[first_index:end,first_index:end]), zdmₚ)::T
                    zi += 1
                end
            end
            if any(tap)
                anneal!(apr, aTeq, dt, tsteps, Tsteps, admₚ)
                ai = 1
                for i ∈ findall(tap)
                    first_index = 1 + floor(Int64,(tinit - crystAge[i])/dt)
                    calcHeAgesₚ[i] = HeAgeSpherical(apatites[ai], @views(Tsteps[first_index:end]), @views(apr[first_index:end,first_index:end]), admₚ)::T
                    ai += 1
                end
            end

            # Calculate log likelihood of proposal
            llnaₚ = diff_ll(Tsteps, dTmax, dTmax_sigma)
            llnaₚ += loglikelihood(admₚ, adm₀) + loglikelihood(zdmₚ, zdm₀)
            simplified && (llnaₚ += -log(npointsₚ))
            σₐ .= simannealsigma.(n, HeAge_sigma, σmodel, σannealing, λannealing)
            llₚ = normpdf_ll(HeAge, σₐ, calcHeAgesₚ) + llnaₚ
            llₗ = normpdf_ll(HeAge, σₐ, calcHeAges) + llna # Recalulate last one too with new σₐ

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - llₗ)

                # Update jumping distribution based on size of current accepted p_move
                if dynamicjumping && r < p_move
                    if agepointsₚ[k] != agepoints[k]
                        σⱼt[k] = ℯ * abs(agepointsₚ[k] - agepoints[k])
                    end
                    if Tpointsₚ[k] != Tpoints[k]
                        σⱼT[k] = ℯ * abs(Tpointsₚ[k] - Tpoints[k])
                    end
                end

                # Update the currently accepted proposal
                adm = admₚ
                zdm = zdmₚ
                ll = llₚ
                llna = llnaₚ
                npoints = npointsₚ
                copyto!(agepoints, 1, agepointsₚ, 1, npoints)
                copyto!(Tpoints, 1, Tpointsₚ, 1, npoints)
                copyto!(constraint.agepoints, constraint.agepointsₚ)
                copyto!(constraint.Tpoints, constraint.Tpointsₚ)
                copyto!(boundary.Tpoints, boundary.Tpointsₚ)
                copyto!(calcHeAges, calcHeAgesₚ)

                # Not critical to the function of the MCMC loop, but critical for recording stationary distribution!
                acceptancedist[n] = true
            end

            # Record results for analysis and troubleshooting
            lldist[n] = llna + normpdf_ll(HeAge, σ, calcHeAges) # Recalculated to constant baseline
            ndist[n] = npoints # distribution of # of points
            σⱼtdist[n] = σⱼt[k]
            σⱼTdist[n] = σⱼT[k]
            HeAgedist[:,n] .= calcHeAges # distribution of He ages

            # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
            collectto!(view(tpointdist, :, n), view(agepoints, 1:npoints), boundary.agepoints, constraint.agepoints)
            collectto!(view(Tpointdist, :, n), view(Tpoints, 1:npoints), boundary.Tpoints, constraint.Tpoints)

            # Update progress meter every `progress_interval` steps
            (mod(n, progress_interval) == 0) && update!(progress, n)
        end
        return (tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist, σⱼtdist, σⱼTdist)
    end
    export MCMC_varkinetics

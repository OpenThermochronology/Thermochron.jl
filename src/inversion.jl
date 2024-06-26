
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


    """
    ```julia
    MCMC(data::NamedTuple, model::NamedTuple, npoints::Int, agepoints::Vector, Tpoints::Vector, unconf::Unconformity, boundary::Boundary, [detail::DetailInterval])
    ```
    Markov chain Monte Carlo time-Temperature inversion of the data specified in `data` and model parameters specified by `model`.

    Returns a tuple of distributions `(tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist)`

    ## Examples
    ```julia
    tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist = MCMC(data, model, npoints, agepoints, Tpoints, unconf, boundary)
    ```
    """
    function MCMC(data::NamedTuple, model::NamedTuple, npoints::Int, agepoints::DenseVector{T}, Tpoints::DenseVector{T}, boundary::Boundary{T}, unconf::Unconformity{T}) where T <: AbstractFloat
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
        Sm = (haskey(data, :Sm) ? T.(data.Th) : zeros(T, size(U)))::DenseVector{T}
        nsteps = model.nsteps::Int
        maxpoints = model.maxpoints::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        totalpoints = maxpoints + boundary.npoints + unconf.npoints::Int
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10.)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : 5.)::T
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

        # Damage models for each mineral
        zdm = ZRDAAM()
        # zdm = ZRDAAM(
        #     DzEa=T(model.DzEa), 
        #     DzD0=T(model.DzD0), 
        #     DN17Ea=T(model.DN17Ea), 
        #     DN17D0=T(model.DN17D0)
        # )
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
        σⱼt = tinit/60
        σⱼT = Tinit/60
        k = 1 # Index of chosen t-T point

        # Proposal probabilities (must sum to 1)
        move = 0.64
        birth = 0.15 # Must equal death
        death = 0.15 # Must equal birth
        movebounds = 0.06

        progress = Progress(nsteps, dt=1, desc="Running MCMC ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps
            @label restart

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
                @label move

                # Move the age of one model point
                agepointsₚ[k] += randn() * σⱼt
                if agepointsₚ[k] < (0 + dt)
                    # Reflecting boundary condition at (0 + dt)
                    agepointsₚ[k] = (0 + dt) - (agepointsₚ[k] - (0 + dt))
                elseif agepointsₚ[k] > (tinit - dt)
                    # Reflecting boundary condition at (tinit - dt)
                    agepointsₚ[k] = (tinit - dt) - (agepointsₚ[k] - (tinit - dt))
                end

                # Move the Temperature of one model point
                Tpointsₚ[k] += randn() * σⱼT
                if Tpointsₚ[k] < Tnow
                    # Reflecting boundary conditions at Tnow (0)
                    Tpointsₚ[k] = Tnow - (Tpointsₚ[k] - Tnow)
                elseif Tpointsₚ[k] > Tinit
                    # Reflecting boundary conditions at Tinit
                    Tpointsₚ[k] = Tinit - (Tpointsₚ[k] - Tinit)
                end

                # Try again if reflected out-of-bounds
                (dt < agepointsₚ[k] < tinit-dt) || @goto move
                (Tnow < Tpointsₚ[k] < Tinit) || @goto move

            elseif (r < move+birth) && (npoints < maxpoints)
                # Birth: add a new model point
                npointsₚ = npoints + 1
                agepointsₚ[npointsₚ] = 0 + rand()*(tinit-0)
                Tpointsₚ[npointsₚ] = Tnow + rand()*(Tinit-Tnow)

            elseif (r < move+birth+death) && (r >= move+birth) && (npoints > minpoints)
                # Death: remove a model point
                npointsₚ = npoints - 1 # Delete last point in array from proposal
                k = ceil(Int, rand()*npoints) # Choose point to delete
                agepointsₚ[k] = agepointsₚ[npoints]
                Tpointsₚ[k] = Tpointsₚ[npoints]

            else
                # Move the temperatures of the starting and ending boundaries
                @. boundary.Tpointsₚ = boundary.T₀ + rand()*boundary.ΔT

                # If there's an imposed unconformity, adjust within parameters
                if unconf.npoints > 0
                    @. unconf.agepointsₚ = unconf.Age₀ + rand()*unconf.ΔAge
                    @. unconf.Tpointsₚ = unconf.T₀ + rand()*unconf.ΔT
                end
            end

            # Recalculate interpolated proposed t-T path
            ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, unconf.agepointsₚ)::StridedVector{T}
            temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, unconf.Tpointsₚ)::StridedVector{T}
            linterp1s!(Tsteps, knot_index, ages, temperatures, agesteps)

            # (maxdiff(Tsteps) > dTmax) && @goto restart
               
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

                # Update jumping distribution based on size of current accepted move
                if dynamicjumping && r < move
                    if agepointsₚ[k] != agepoints[k]
                        σⱼt = max(ℯ * abs(agepointsₚ[k] - agepoints[k]), dt)
                    end
                    if Tpointsₚ[k] != Tpoints[k]
                        σⱼT = max(ℯ * abs(Tpointsₚ[k] - Tpoints[k]), (Tinit-Tnow)/100)
                    end
                end

                # Update the currently accepted proposal
                ll = llₚ
                llna = llnaₚ
                npoints = npointsₚ
                copyto!(agepoints, 1, agepointsₚ, 1, npoints)
                copyto!(Tpoints, 1, Tpointsₚ, 1, npoints)
                copyto!(unconf.agepoints, unconf.agepointsₚ)
                copyto!(unconf.Tpoints, unconf.Tpointsₚ)
                copyto!(boundary.Tpoints, boundary.Tpointsₚ)
                copyto!(calcHeAges, calcHeAgesₚ)

                # Not critical to the function of the MCMC loop, but critical for recording stationary distribution!
                acceptancedist[n] = true
            end

            # Record results for analysis and troubleshooting
            lldist[n] = llna + normpdf_ll(HeAge, σ, calcHeAges) # Recalculated to constant baseline
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
    function MCMC(data::NamedTuple, model::NamedTuple, npoints::Int, agepoints::DenseVector{T}, Tpoints::DenseVector{T}, boundary::Boundary{T}, unconf::Unconformity{T}, detail::DetailInterval{T}) where T <: AbstractFloat
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
        Sm = (haskey(data, :Sm) ? T.(data.Th) : zeros(T, size(U)))::DenseVector{T}
        nsteps = model.nsteps::Int
        maxpoints = model.maxpoints::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        totalpoints = maxpoints + boundary.npoints + unconf.npoints::Int
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10.)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : 5.)::T
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
        σⱼt = tinit/60
        σⱼT = Tinit/60
        k = 1 # Index of chosen t-T point

        # Proposal probabilities (must sum to 1)
        move = 0.64
        birth = 0.15 # Must equal death
        death = 0.15 # Must equal birth
        movebounds = 0.06

        progress = Progress(nsteps, dt=1, desc="Running MCMC ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps
            enoughpoints = min(pointsininterval(agepoints, npoints, detail.agemin, detail.agemax, dt), detail.minpoints)::Int
            @label restart

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
                @label move

                # Move the age of one model point
                agepointsₚ[k] += randn() * σⱼt
                if agepointsₚ[k] < (0 + dt)
                    # Reflecting boundary condition at (0 + dt)
                    agepointsₚ[k] = (0 + dt) - (agepointsₚ[k] - (0 + dt))
                elseif agepointsₚ[k] > (tinit - dt)
                    # Reflecting boundary condition at (tinit - dt)
                    agepointsₚ[k] = (tinit - dt) - (agepointsₚ[k] - (tinit - dt))
                end

                # Move the Temperature of one model point
                Tpointsₚ[k] += randn() * σⱼT
                if Tpointsₚ[k] < Tnow
                    # Reflecting boundary conditions at Tnow (0)
                    Tpointsₚ[k] = Tnow - (Tpointsₚ[k] - Tnow)
                elseif Tpointsₚ[k] > Tinit
                    # Reflecting boundary conditions at Tinit
                    Tpointsₚ[k] = Tinit - (Tpointsₚ[k] - Tinit)
                end

                # Try again if reflected out-of-bounds
                (dt < agepointsₚ[k] < tinit-dt) || @goto move
                (Tnow < Tpointsₚ[k] < Tinit) || @goto move

            elseif (r < move+birth) && (npoints < maxpoints)
                # Birth: add a new model point
                npointsₚ = npoints + 1
                agepointsₚ[npointsₚ] = 0 + rand()*(tinit-0)
                Tpointsₚ[npointsₚ] = Tnow + rand()*(Tinit-Tnow)

            elseif (r < move+birth+death) && (r >= move+birth) && (npoints > max(minpoints, detail.minpoints))
                # Death: remove a model point
                npointsₚ = npoints - 1 # Delete last point in array from proposal
                k = ceil(Int, rand()*npoints) # Choose point to delete
                agepointsₚ[k] = agepointsₚ[npoints]
                Tpointsₚ[k] = Tpointsₚ[npoints]

            else
                # Move the temperatures of the starting and ending boundaries
                @. boundary.Tpointsₚ = boundary.T₀ + rand()*boundary.ΔT

                # If there's an imposed unconformity, adjust within parameters
                if unconf.npoints > 0
                    @. unconf.agepointsₚ = unconf.Age₀ + rand()*unconf.ΔAge
                    @. unconf.Tpointsₚ = unconf.T₀ + rand()*unconf.ΔT
                end
            end

            # Recalculate interpolated proposed t-T path
            ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, unconf.agepointsₚ)::StridedVector{T}
            temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, unconf.Tpointsₚ)::StridedVector{T}
            linterp1s!(Tsteps, knot_index, ages, temperatures, agesteps)
        
            # (maxdiff(Tsteps) > dTmax) && @goto restart
            (pointsininterval(agepointsₚ, npointsₚ, detail.agemin, detail.agemax, dt) < enoughpoints) && @goto restart
                    
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

                # Update jumping distribution based on size of current accepted move
                if dynamicjumping && r < move
                    if agepointsₚ[k] != agepoints[k]
                        σⱼt = max(ℯ * abs(agepointsₚ[k] - agepoints[k]), dt)
                    end
                    if Tpointsₚ[k] != Tpoints[k]
                        σⱼT = max(ℯ * abs(Tpointsₚ[k] - Tpoints[k]), (Tinit-Tnow)/100)
                    end
                end

                # Update the currently accepted proposal
                ll = llₚ
                llna = llnaₚ
                npoints = npointsₚ
                copyto!(agepoints, 1, agepointsₚ, 1, npoints)
                copyto!(Tpoints, 1, Tpointsₚ, 1, npoints)
                copyto!(unconf.agepoints, unconf.agepointsₚ)
                copyto!(unconf.Tpoints, unconf.Tpointsₚ)
                copyto!(boundary.Tpoints, boundary.Tpointsₚ)
                copyto!(calcHeAges, calcHeAgesₚ)

                # Not critical to the function of the MCMC loop, but critical for recording stationary distribution!
                acceptancedist[n] = true
            end

            # Record results for analysis and troubleshooting
            lldist[n] = llna + normpdf_ll(HeAge, σ, calcHeAges) # Recalculated to constant baseline
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
    export MCMC

    function MCMC_varkinetics(data::NamedTuple, model::NamedTuple, npoints::Int, agepoints::DenseVector{T}, Tpoints::DenseVector{T}, boundary::Boundary{T}, unconf::Unconformity{T}, detail::DetailInterval{T}) where T <: AbstractFloat
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
        Sm = (haskey(data, :Sm) ? T.(data.Th) : zeros(T, size(U)))::DenseVector{T}
        nsteps = model.nsteps::Int
        maxpoints = model.maxpoints::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        totalpoints = maxpoints + boundary.npoints + unconf.npoints::Int
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10.)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : 5.)::T
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
        # for kinetics
        σaD0 = adm.D0L_logsigma
        σaEa = adm.EaL_logsigma
        σzD0 = (zdm.DzD0_logsigma+zdm.DN17D0_logsigma)/2
        σzEa = (zdm.DzEa_logsigma+zdm.DN17Ea_logsigma)/2
        # for temperature and time
        σⱼt = tinit/60
        σⱼT = Tinit/60
        k = 1 # Index of chosen t-T point

        # Proposal probabilities (must sum to 1)
        move = 0.6
        birth = 0.15 # Must equal death
        death = 0.15 # Must equal birth
        movebounds = 0.06
        movekinetics = 0.04

        progress = Progress(nsteps, dt=1, desc="Running MCMC ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps
            enoughpoints = min(pointsininterval(agepoints, npoints, detail.agemin, detail.agemax, dt), detail.minpoints)::Int
            @label restart

            # Copy proposal from last accepted solution
            admₚ = adm
            zdmₚ = zdm
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
                @label move

                # Move the age of one model point
                agepointsₚ[k] += randn() * σⱼt
                if agepointsₚ[k] < (0 + dt)
                    # Reflecting boundary condition at (0 + dt)
                    agepointsₚ[k] = (0 + dt) - (agepointsₚ[k] - (0 + dt))
                elseif agepointsₚ[k] > (tinit - dt)
                    # Reflecting boundary condition at (tinit - dt)
                    agepointsₚ[k] = (tinit - dt) - (agepointsₚ[k] - (tinit - dt))
                end

                # Move the Temperature of one model point
                Tpointsₚ[k] += randn() * σⱼT
                if Tpointsₚ[k] < Tnow
                    # Reflecting boundary conditions at Tnow (0)
                    Tpointsₚ[k] = Tnow - (Tpointsₚ[k] - Tnow)
                elseif Tpointsₚ[k] > Tinit
                    # Reflecting boundary conditions at Tinit
                    Tpointsₚ[k] = Tinit - (Tpointsₚ[k] - Tinit)
                end

                # Try again if reflected out-of-bounds
                (dt < agepointsₚ[k] < tinit-dt) || @goto move
                (Tnow < Tpointsₚ[k] < Tinit) || @goto move

            elseif (r < move+birth) && (npoints < maxpoints)
                # Birth: add a new model point
                npointsₚ = npoints + 1
                agepointsₚ[npointsₚ] = 0 + rand()*(tinit-0)
                Tpointsₚ[npointsₚ] = Tnow + rand()*(Tinit-Tnow)

            elseif (r < move+birth+death) && (r >= move+birth) && (npoints > max(minpoints, detail.minpoints))
                # Death: remove a model point
                npointsₚ = npoints - 1 # Delete last point in array from proposal
                k = ceil(Int, rand()*npoints) # Choose point to delete
                agepointsₚ[k] = agepointsₚ[npoints]
                Tpointsₚ[k] = Tpointsₚ[npoints]

            elseif (r < move+birth+death+movebounds)
                # Move the temperatures of the starting and ending boundaries
                @. boundary.Tpointsₚ = boundary.T₀ + rand()*boundary.ΔT

                # If there's an imposed unconformity, adjust within parameters
                if unconf.npoints > 0
                    @. unconf.agepointsₚ = unconf.Age₀ + rand()*unconf.ΔAge
                    @. unconf.Tpointsₚ = unconf.T₀ + rand()*unconf.ΔT
                end

            else
                # Adjust kinetic parameters, one at a time
                if any(tzr)
                    rn = rand(1:4)
                    zdmₚ = ZRDAAM(
                        DzEa = (rn==1) ? exp(log(zdm.DzEa)+randn()*σzEa) : zdm.DzEa,
                        DzD0 = (rn==2) ? exp(log(zdm.DzD0)+randn()*σzD0) : zdm.DzD0,
                        DN17Ea = (rn==3) ? exp(log(zdm.DN17Ea)+randn()*σzEa) : zdm.DN17Ea,
                        DN17D0 = (rn==4) ? exp(log(zdm.DN17D0)+randn()*σzD0) : zdm.DN17D0,
                    )
                end
                if any(tap)
                    rn = rand(1:3)
                    admₚ = RDAAM(
                        D0L = (rn==1) ? exp(log(adm.D0L)+randn()*σaD0) : adm.D0L,
                        EaL = (rn==2) ? exp(log(adm.EaL)+randn()*σaEa) : adm.EaL,
                        EaTrap = (rn==3) ? exp(log(adm.EaTrap)+randn()*σaEa) : adm.EaTrap,
                    )
                end
            end

            # Recalculate interpolated proposed t-T path
            ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, unconf.agepointsₚ)::StridedVector{T}
            temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, unconf.Tpointsₚ)::StridedVector{T}
            linterp1s!(Tsteps, knot_index, ages, temperatures, agesteps)
    
            # (maxdiff(Tsteps) > dTmax) && @goto restart
            (pointsininterval(agepointsₚ, npointsₚ, detail.agemin, detail.agemax, dt) < enoughpoints) && @goto restart
               
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

                # Update jumping distribution based on size of current accepted move
                if dynamicjumping && r < move
                    if agepointsₚ[k] != agepoints[k]
                        σⱼt = max(ℯ * abs(agepointsₚ[k] - agepoints[k]), dt)
                    end
                    if Tpointsₚ[k] != Tpoints[k]
                        σⱼT = max(ℯ * abs(Tpointsₚ[k] - Tpoints[k]), (Tinit-Tnow)/100)
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
                copyto!(unconf.agepoints, unconf.agepointsₚ)
                copyto!(unconf.Tpoints, unconf.Tpointsₚ)
                copyto!(boundary.Tpoints, boundary.Tpointsₚ)
                copyto!(calcHeAges, calcHeAgesₚ)

                # Not critical to the function of the MCMC loop, but critical for recording stationary distribution!
                acceptancedist[n] = true
            end

            # Record results for analysis and troubleshooting
            lldist[n] = llna + normpdf_ll(HeAge, σ, calcHeAges) # Recalculated to constant baseline
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
    function MCMC_varkinetics(data::NamedTuple, model::NamedTuple, npoints::Int, agepoints::DenseVector{T}, Tpoints::DenseVector{T}, boundary::Boundary{T}, unconf::Unconformity{T}) where T <: AbstractFloat
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
        Sm = (haskey(data, :Sm) ? T.(data.Th) : zeros(T, size(U)))::DenseVector{T}
        nsteps = model.nsteps::Int
        maxpoints = model.maxpoints::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        totalpoints = maxpoints + boundary.npoints + unconf.npoints::Int
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10.)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : 5.)::T
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
        # for kinetics
        σaD0 = adm.D0L_logsigma
        σaEa = adm.EaL_logsigma
        σzD0 = (zdm.DzD0_logsigma+zdm.DN17D0_logsigma)/2
        σzEa = (zdm.DzEa_logsigma+zdm.DN17Ea_logsigma)/2
        # for temperature and time
        σⱼt = tinit/60
        σⱼT = Tinit/60
        k = 1 # Index of chosen t-T point
    
        # Proposal probabilities (must sum to 1)
        move = 0.6
        birth = 0.15 # Must equal death
        death = 0.15 # Must equal birth
        movebounds = 0.06
        movekinetics = 0.04
    
        progress = Progress(nsteps, dt=1, desc="Running MCMC ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps
            @label restart
    
            # Copy proposal from last accepted solution
            admₚ = adm
            zdmₚ = zdm
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
                @label move
    
                # Move the age of one model point
                agepointsₚ[k] += randn() * σⱼt
                if agepointsₚ[k] < (0 + dt)
                    # Reflecting boundary condition at (0 + dt)
                    agepointsₚ[k] = (0 + dt) - (agepointsₚ[k] - (0 + dt))
                elseif agepointsₚ[k] > (tinit - dt)
                    # Reflecting boundary condition at (tinit - dt)
                    agepointsₚ[k] = (tinit - dt) - (agepointsₚ[k] - (tinit - dt))
                end
    
                # Move the Temperature of one model point
                Tpointsₚ[k] += randn() * σⱼT
                if Tpointsₚ[k] < Tnow
                    # Reflecting boundary conditions at Tnow (0)
                    Tpointsₚ[k] = Tnow - (Tpointsₚ[k] - Tnow)
                elseif Tpointsₚ[k] > Tinit
                    # Reflecting boundary conditions at Tinit
                    Tpointsₚ[k] = Tinit - (Tpointsₚ[k] - Tinit)
                end
    
                # Try again if reflected out-of-bounds
                (dt < agepointsₚ[k] < tinit-dt) || @goto move
                (Tnow < Tpointsₚ[k] < Tinit) || @goto move
                
            elseif (r < move+birth) && (npoints < maxpoints)
                # Birth: add a new model point
                npointsₚ = npoints + 1
                agepointsₚ[npointsₚ] = 0 + rand()*(tinit-0)
                Tpointsₚ[npointsₚ] = Tnow + rand()*(Tinit-Tnow)
    
            elseif (r < move+birth+death) && (r >= move+birth) && (npoints > minpoints)
                # Death: remove a model point
                npointsₚ = npoints - 1 # Delete last point in array from proposal
                k = ceil(Int, rand()*npoints) # Choose point to delete
                agepointsₚ[k] = agepointsₚ[npoints]
                Tpointsₚ[k] = Tpointsₚ[npoints]
    
            elseif (r < move+birth+death+movebounds)
                # Move the temperatures of the starting and ending boundaries
                @. boundary.Tpointsₚ = boundary.T₀ + rand()*boundary.ΔT

                # If there's an imposed unconformity, adjust within parameters
                if unconf.npoints > 0
                    @. unconf.agepointsₚ = unconf.Age₀ + rand()*unconf.ΔAge
                    @. unconf.Tpointsₚ = unconf.T₀ + rand()*unconf.ΔT
                end
    
            else
                # Adjust kinetic parameters, one at a time
                if any(tzr)
                    rn = rand(1:4)
                    zdmₚ = ZRDAAM(
                        DzEa = (rn==1) ? exp(log(zdm.DzEa)+randn()*σzEa) : zdm.DzEa,
                        DzD0 = (rn==2) ? exp(log(zdm.DzD0)+randn()*σzD0) : zdm.DzD0,
                        DN17Ea = (rn==3) ? exp(log(zdm.DN17Ea)+randn()*σzEa) : zdm.DN17Ea,
                        DN17D0 = (rn==4) ? exp(log(zdm.DN17D0)+randn()*σzD0) : zdm.DN17D0,
                    )
                end
                if any(tap)
                    rn = rand(1:3)
                    admₚ = RDAAM(
                        D0L = (rn==1) ? exp(log(adm.D0L)+randn()*σaD0) : adm.D0L,
                        EaL = (rn==2) ? exp(log(adm.EaL)+randn()*σaEa) : adm.EaL,
                        EaTrap = (rn==3) ? exp(log(adm.EaTrap)+randn()*σaEa) : adm.EaTrap,
                    )
                end
            end

            # Recalculate interpolated proposed t-T path
            ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, unconf.agepointsₚ)::StridedVector{T}
            temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, unconf.Tpointsₚ)::StridedVector{T}
            linterp1s!(Tsteps, knot_index, ages, temperatures, agesteps)
        
            # (maxdiff(Tsteps) > dTmax) && @goto restart

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
    
                # Update jumping distribution based on size of current accepted move
                if dynamicjumping && r < move
                    if agepointsₚ[k] != agepoints[k]
                        σⱼt = max(ℯ * abs(agepointsₚ[k] - agepoints[k]), dt)
                    end
                    if Tpointsₚ[k] != Tpoints[k]
                        σⱼT = max(ℯ * abs(Tpointsₚ[k] - Tpoints[k]), (Tinit-Tnow)/100)
                    end
                end
    
                # Update the currently accepted proposal
                adm = admₚ
                zdm = zdmₚ
                ll = llₚ
                llna = llnaₚ
                npoints = npointsₚ
                copyto!(agepoints, agepointsₚ)
                copyto!(Tpoints, Tpointsₚ)
                copyto!(unconf.agepoints, unconf.agepointsₚ)
                copyto!(unconf.Tpoints, unconf.Tpointsₚ)
                copyto!(boundary.Tpoints, boundary.Tpointsₚ)
                copyto!(calcHeAges, calcHeAgesₚ)
    
                # Not critical to the function of the MCMC loop, but critical for recording stationary distribution!
                acceptancedist[n] = true
            end
    
            # Record results for analysis and troubleshooting
            lldist[n] = llna + normpdf_ll(HeAge, σ, calcHeAges) # Recalculated to constant baseline
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
    export MCMC_varkinetics

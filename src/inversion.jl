
    """
    ```julia
    simannealsigma(n, sigma_analytical; [simannealmodel::NamedTuple])
    ```
    To avoid getting stuck in local optima, decrease uncertainty slowly by
    simulated annealing. Parameters are specified as a tuple `simannealmodel` of the
    form (σₘ, σᵢ, λ), where annealing uncertainty declines from `σᵢ+σₘ` to `σₘ`
    with a decay constant of λ.

    Returns the annealing uncertainty added in quadrature with analytical
    uncertainty, or in other words

        sigma_annealing = σᵢ*exp(-λ*n) + σₘ
        sigma = sqrt(sigma_analytical^2 + sigma_annealing^2)

    """
    function simannealsigma(n::Number, sigma_analytical::Number; simannealmodel::NamedTuple=(σModel=25.0, σAnnealing=35.0, λAnnealing=10/10^5))
        mdl = simannealmodel
        sigma_combined = mdl.σAnnealing * exp(-mdl.λAnnealing*n) + mdl.σModel
        return sqrt(sigma_analytical^2 + sigma_combined^2)
    end
    export simannealsigma

    # Utility function for agePoint and TPoint buffers
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

    # Utility function for checking maximum reheating rate
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

    """
    ```julia
    MCMC_vartcryst(data::NamedTuple, model::NamedTuple, nPoints::Int, agePoints::Vector, TPoints::Vector, unconf::NamedTuple, boundary::NamedTuple)
    ```
    Markov chain Monte Carlo time-Temperature inversion of the data specified in `data` and model parameters specified by `model`.

    Returns a tuple of distributions `(TStepdist, HeAgedist, ndist, lldist, acceptancedist)`

    ## Examples
    ```julia
    TStepdist, HeAgedist, ndist, lldist, acceptancedist = MCMC_vartcryst(data, model, nPoints, agePoints, TPoints, unconf, boundary)
    ```
    """
    function MCMC_vartcryst(data::NamedTuple, model::NamedTuple, nPoints::Int, agePoints::AbstractVector, TPoints::AbstractVector, unconf::NamedTuple, boundary::NamedTuple)
        @assert firstindex(agePoints) === 1
        @assert firstindex(TPoints) === 1
        nsteps = model.nsteps::Int
        maxPoints = model.maxPoints::Int
        minPoints = (haskey(model, :minPoints) ? model.minPoints : 1)::Int

        # Calculate number of boundary and unconformity points and allocate buffer for interpolating
        extraPoints = length(boundary.agePoints) + length(unconf.agePoints)
        agePointBuffer = similar(agePoints, maxPoints+extraPoints)
        TPointBuffer = similar(agePoints, maxPoints+extraPoints)
        knot_index = similar(model.ageSteps, Int)

        # Calculate model ages for initial proposal
        ages = collectto!(agePointBuffer, view(agePoints, 1:nPoints), boundary.agePoints, unconf.agePoints)
        temperatures = collectto!(TPointBuffer, view(TPoints, 1:nPoints), boundary.TPoints, unconf.TPoints)
        TSteps = linterp1s(ages, temperatures, model.ageSteps)
        calcHeAges = Array{Float64}(undef, size(data.HeAge))
        pr, Teq = anneal(model.dt, model.tSteps, TSteps, ZRDAAM()) # Damage annealing history

        # Prepare a Mineral object for each analysis
        diffusionmodel = (DzEa=model.DzEa, DzD0=model.DzD0, DN17Ea=model.DN17Ea, DN17D0=model.DN17D0)
        zircons = Array{Zircon{Float64}}(undef, length(data.halfwidth))
        for i=1:length(zircons)
            # Iterate through each grain, calculate the modeled age for each
            first_index = 1 + floor(Int64,(model.tInit - data.CrystAge[i])/model.dt)
            zircons[i] = Zircon(data.halfwidth[i], model.dr, data.U[i], data.Th[i], model.dt, model.ageSteps[first_index:end])
            calcHeAges[i] = HeAgeSpherical(zircons[i], @views(TSteps[first_index:end]), @views(pr[first_index:end,first_index:end]), diffusionmodel)
        end

        # Simulated annealing of uncertainty
        simannealmodel = (σModel=model.σModel, σAnnealing=model.σAnnealing, λAnnealing=model.λAnnealing)
        σₐ = simannealsigma.(1, data.HeAge_sigma; simannealmodel)
        σₙ = simannealsigma.(nsteps, data.HeAge_sigma; simannealmodel)

        # Log-likelihood for initial proposal
        ll = normpdf_ll(data.HeAge, σₐ, calcHeAges)
        if model.simplified
            ll -= log(nPoints)
        end

        # Variables to hold proposals
        llₚ = ll
        nPointsₚ = nPoints
        agePointsₚ = similar(agePoints)
        TPointsₚ = similar(TPoints)
        calcHeAgesₚ = similar(calcHeAges)

        # distributions to populate
        HeAgedist = Array{Float64}(undef, length(data.HeAge), nsteps)
        TStepdist = Array{Float64}(undef, length(model.tSteps), nsteps)
        lldist = Array{Float64}(undef, nsteps)
        ndist = zeros(Int, nsteps)
        acceptancedist = zeros(Bool, nsteps)

        # Standard deviations of Gaussian proposal distributions for temperature and time
        t_sigma = model.tInit/60
        T_sigma = model.TInit/60

        # Proposal probabilities (must sum to 1)
        move = 0.64
        birth = 0.15
        death = 0.15 # Should equal birth
        movebounds = 0.06

        progress = Progress(nsteps, dt=1, desc="Running MCMC ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps

            # Copy proposal from last accepted solution
            nPointsₚ = nPoints
            copyto!(agePointsₚ, agePoints)
            copyto!(TPointsₚ, TPoints)
            copyto!(unconf.agePointsₚ, unconf.agePoints)
            copyto!(unconf.TPointsₚ, unconf.TPoints)
            copyto!(boundary.TPointsₚ, boundary.TPoints)

            # Adjust the proposal
            r = rand()
            if r < move
                # Move one t-T point
                k = ceil(Int, rand() * nPoints)
                while true

                # Move the age of one model point
                agePointsₚ[k] += randn() * t_sigma
                if agePointsₚ[k] < (0 + model.dt)
                    # Don't let any point get too close to 0
                    agePointsₚ[k] = 0 + model.dt
                elseif agePointsₚ[k] > (model.tInit - model.dt)
                    # Don't let any point get too close to tInit
                    agePointsₚ[k] = model.tInit - model.dt
                end

                # Move the Temperature of one model point
                TPointsₚ[k] += randn() * T_sigma
                if TPointsₚ[k] < 0
                    # Don't allow T<0
                    TPointsₚ[k] = 0
                elseif TPointsₚ[k] > model.TInit
                    # Don't allow T>TInit
                    TPointsₚ[k] = model.TInit
                end

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agePointBuffer, view(agePointsₚ, 1:nPointsₚ), boundary.agePoints, unconf.agePointsₚ)
                temperatures = collectto!(TPointBuffer, view(TPointsₚ, 1:nPointsₚ), boundary.TPointsₚ, unconf.TPointsₚ)
                linterp1s!(TSteps, knot_index, ages, temperatures, model.ageSteps)

                # Retry unless we have satisfied the maximum reheating rate
                (maxdiff(TSteps) < model.dTmax) && break
                end
            elseif (r < move+birth) && (nPointsₚ < maxPoints)
                # Birth: add a new model point
                nPointsₚ += 1
                while true
                agePointsₚ[nPointsₚ] = rand()*model.tInit
                TPointsₚ[nPointsₚ] = model.TNow + rand()*(model.TInit-model.TNow)

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agePointBuffer, view(agePointsₚ, 1:nPointsₚ), boundary.agePoints, unconf.agePointsₚ)
                temperatures = collectto!(TPointBuffer, view(TPointsₚ, 1:nPointsₚ), boundary.TPointsₚ, unconf.TPointsₚ)
                linterp1s!(TSteps, knot_index, ages, temperatures, model.ageSteps)

                # Retry unless we have satisfied the maximum reheating rate
                (maxdiff(TSteps) < model.dTmax) && break
                end
            elseif (r < move+birth+death) && (r >= move+birth) && (nPointsₚ > minPoints)
                # Death: remove a model point
                nPointsₚ -= 1 # Delete last point in array from proposal
                while true
                k = ceil(Int, rand()*nPoints) # Choose point to delete
                agePointsₚ[k] = agePointsₚ[nPoints]
                TPointsₚ[k] = TPointsₚ[nPoints]

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agePointBuffer, view(agePointsₚ, 1:nPointsₚ), boundary.agePoints, unconf.agePointsₚ)
                temperatures = collectto!(TPointBuffer, view(TPointsₚ, 1:nPointsₚ), boundary.TPointsₚ, unconf.TPointsₚ)
                linterp1s!(TSteps, knot_index, ages, temperatures, model.ageSteps)

                # Retry unless we have satisfied the maximum reheating rate
                (maxdiff(TSteps) < model.dTmax) && break
                end
            else
                # Move boundary conditions
                while true
                # Move the temperatures of the starting and ending boundaries
                @. boundary.TPointsₚ = boundary.T₀ + rand()*boundary.ΔT

                # If there's an imposed unconformity, adjust within parameters
                if length(unconf.agePoints) > 0
                    @. unconf.agePointsₚ = unconf.Age₀ + rand()*unconf.ΔAge
                    @. unconf.TPointsₚ = unconf.T₀ + rand()*unconf.ΔT
                end

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agePointBuffer, view(agePointsₚ, 1:nPointsₚ), boundary.agePoints, unconf.agePointsₚ)
                temperatures = collectto!(TPointBuffer, view(TPointsₚ, 1:nPointsₚ), boundary.TPointsₚ, unconf.TPointsₚ)
                linterp1s!(TSteps, knot_index, ages, temperatures, model.ageSteps)

                # Retry unless we have satisfied the maximum reheating rate
                (maxdiff(TSteps) < model.dTmax) && break
                end
            end

             # Calculate model ages for each grain
            anneal!(pr, Teq, model.dt, model.tSteps, TSteps, ZRDAAM())
            for i=1:length(zircons)
                first_index = 1 + floor(Int64,(model.tInit - data.CrystAge[i])/model.dt)
                calcHeAgesₚ[i] = HeAgeSpherical(zircons[i], @views(TSteps[first_index:end]), @views(pr[first_index:end,first_index:end]), diffusionmodel)
            end

            # Calculate log likelihood of proposal
            σₐ .= simannealsigma.(n, data.HeAge_sigma; simannealmodel)
            llₚ = normpdf_ll(data.HeAge, σₐ, calcHeAgesₚ)
            llₗ = normpdf_ll(data.HeAge, σₐ, calcHeAges) # Recalulate last one too with new σₐ
            if model.simplified # slightly penalize more complex t-T paths
                llₚ -= log(nPointsₚ)
                llₗ -= log(nPoints)
            end

            # Accept or reject proposal based on likelihood
            # To avoid numerical problems with diffusion code, also reject proposal
            # if maximum proposed heating rate is greater than 25C per timestep.
            # (Fast cooling should not be a problem, however)
            if log(rand()) < (llₚ - llₗ)
                ll = llₚ
                nPoints = nPointsₚ
                copyto!(agePoints, agePointsₚ)
                copyto!(TPoints, TPointsₚ)
                copyto!(unconf.agePoints, unconf.agePointsₚ)
                copyto!(unconf.TPoints, unconf.TPointsₚ)
                copyto!(boundary.TPoints, boundary.TPointsₚ)
                copyto!(calcHeAges, calcHeAgesₚ)

                # This is saved for ouput, but not critical to the function of the MCMC loop
                acceptancedist[n] = true
            end

            # Record results for analysis and troubleshooting
            lldist[n] = normpdf_ll(data.HeAge, σₙ, calcHeAges) # Recalculated to constant baseline
            ndist[n] = nPoints # distribution of # of points
            HeAgedist[:,n] .= calcHeAges # distribution of He ages

            # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
            TStepdist[:,n] .= TSteps # distribution of T paths

            # Update progress meter every `progress_interval` steps
            (mod(n, progress_interval) == 0) && update!(progress, n)
        end
        return (TStepdist, HeAgedist, ndist, lldist, acceptancedist)
    end
    export MCMC_vartcryst

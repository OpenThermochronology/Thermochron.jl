
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
        i₀ = 1
        copyto!(buffer, i₀, a, 1, length(a))
        i₀ += length(a)
        copyto!(buffer, i₀, b, 1, length(b))
        i₀ += length(b)
        copyto!(buffer, i₀, c, 1, length(c))
        n = length(a) + length(b) + length(c)
        return n
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

    Returns a tuple of distributions `(TStepDist, HeAgeDist, nDist, llDist, acceptanceDist)`

    ## Examples
    ```julia
    TStepDist, HeAgeDist, nDist, llDist, acceptanceDist = MCMC_vartcryst(data, model, nPoints, agePoints, TPoints, unconf, boundary)
    ```
    """
    function MCMC_vartcryst(data::NamedTuple, model::NamedTuple, nPoints::Int, agePoints::AbstractVector, TPoints::AbstractVector, unconf::NamedTuple, boundary::NamedTuple)
        # Calculate model ages for initial proposal
        TSteps = linterp1s([view(agePoints, 1:nPoints) ; boundary.agePoints ; unconf.agePoints],
                           [view(TPoints, 1:nPoints) ; boundary.TPoints ; unconf.TPoints], model.ageSteps)
        CalcHeAges = Array{Float64}(undef, size(data.HeAge))
        pr, Teq = anneal(model.dt, model.tSteps, TSteps, ZRDAAM()) # Damage annealing history
        zircons = Array{Zircon{Float64}}(undef, length(data.halfwidth))
        for i=1:length(data.halfwidth)
            # Iterate through each grain, calculate the modeled age for each
            first_index = 1 + floor(Int64,(model.tInit - data.CrystAge[i])/model.dt)
            zircons[i] = Zircon(data.halfwidth[i], model.dr, data.U[i], data.Th[i], model.dt, model.ageSteps[first_index:end])
            CalcHeAges[i] = HeAgeSpherical(zircons[i], TSteps[first_index:end], pr[first_index:end,first_index:end], model)
        end
        AnnealedSigma = simannealsigma.(1, data.HeAge_sigma; simannealmodel=model)
        UnAnnealedSigma = simannealsigma.(model.nsteps, data.HeAge_sigma; simannealmodel=model)

        # Log-likelihood for initial proposal
        ll = normpdf_ll(data.HeAge, AnnealedSigma, CalcHeAges)
        if model.simplified
            ll -= log(nPoints)
        end

        # Variables to hold proposals
        llₚ = ll
        nPointsₚ = nPoints
        agePointsₚ = similar(agePoints)
        TPointsₚ = similar(TPoints)
        TStepsₚ = similar(TSteps)
        CalcHeAgesₚ = similar(CalcHeAges)

        # Distributions to populate
        HeAgeDist = Array{Float64}(undef, length(data.HeAge), model.nsteps)
        TStepDist = Array{Float64}(undef, length(model.tSteps), model.nsteps)
        llDist = Array{Float64}(undef, model.nsteps)
        nDist = zeros(Int, model.nsteps)
        acceptanceDist = zeros(Bool, model.nsteps)

        # Standard deviations of Gaussian proposal distributions for temperature and time
        t_sigma = model.tInit/60
        T_sigma = model.TInit/60

        # Proposal probabilities (must sum to 1)
        move = 0.64
        birth = 0.15
        death = 0.15 # Should equal birth
        movebounds = 0.06
        maxattempts = 1000

        @showprogress 10 "Running MCMC..." for n=1:model.nsteps

            # Copy proposal from last accepted solution
            nPointsₚ = nPoints
            copyto!(agePointsₚ, agePoints)
            copyto!(TPointsₚ, TPoints)
            copyto!(TStepsₚ, TSteps)
            copyto!(unconf.agePointsₚ, unconf.agePoints)
            copyto!(unconf.TPointsₚ, unconf.TPoints)
            copyto!(boundary.TPointsₚ, boundary.TPoints)

            # Adjust the proposal
            r = rand()
            if r < move
                # Move the age of one model point
                for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                    k = ceil(Int, rand() * nPoints)

                    agePointsₚ[k] += randn() * t_sigma
                    if agePointsₚ[k] < model.dt
                        # Don't let any point get too close to 0
                        agePointsₚ[k] += (model.dt - agePointsₚ[k])
                    elseif agePointsₚ[k] > (model.tInit - model.dt)
                        # Don't let any point get too close to model.tInit
                        agePointsₚ[k] -= (agePointsₚ[k] - (model.tInit - model.dt))
                    end
                    # Move the Temperature of one model point
                    if TPointsₚ[k] < 0
                        # Don't allow T<0
                        TPointsₚ[k] = 0
                    elseif TPointsₚ[k] > model.TInit
                        # Don't allow T>model.TInit
                        TPointsₚ[k] = model.TInit
                    end

                    # Interpolate proposed t-T path
                    TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                        [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], model.ageSteps)

                    # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                    maximum(diff(TStepsₚ)) < model.dTmax && break

                    # Copy last accepted solution to re-modify if we don't break
                    copyto!(agePointsₚ, agePoints)
                    copyto!(TPointsₚ, TPoints)
                end

            elseif (r < move+birth) && (nPointsₚ < model.maxPoints)
                # Birth: add a new model point
                nPointsₚ += 1
                for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                    agePointsₚ[nPointsₚ] = rand()*model.tInit
                    TPointsₚ[nPointsₚ] = rand()*model.TInit

                    # Interpolate proposed t-T path
                    TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                        [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], model.ageSteps)

                    # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                    maximum(diff(TStepsₚ)) < model.dTmax && break
                end

            elseif (r < move+birth+death) && (r >= move+birth) && (nPointsₚ > 1)
                # Death: remove a model point
                nPointsₚ -= 1 # Delete last point in array from proposal
                for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                    k = ceil(Int, rand()*nPoints) # Choose point to delete
                    agePointsₚ[k] = agePointsₚ[nPoints]
                    TPointsₚ[k] = TPointsₚ[nPoints]

                    # Interpolate proposed t-T path
                    TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                        [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], model.ageSteps)

                    # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                    maximum(diff(TStepsₚ)) < model.dTmax && break
                end

            else
                # Move boundary conditions
                for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                    if rand() < 0.5
                        # Allow the present temperature to vary from 0 to 10 degrees C
                        boundary.TPointsₚ[1] = 0+rand()*10
                    else
                        # Allow the initial temperature to vary from model.TInit to model.TInit-50 C
                        boundary.TPointsₚ[2] = model.TInit-rand()*50
                    end
                    if length(unconf.agePointsₚ) > 0
                        # If there's an imposed unconformity, adjust within parameters
                        @. unconf.agePointsₚ = unconf.Age₀ + rand()*unconf.ΔAge
                        @. unconf.TPointsₚ = unconf.T₀ + rand()*unconf.ΔT
                    end

                    # Recalculate interpolated proposed t-T path
                    TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                            [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], model.ageSteps)

                    # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                    maximum(diff(TStepsₚ)) < model.dTmax && break

                    # Copy last accepted solution to re-modify if we don't break
                    copyto!(unconf.agePointsₚ, unconf.agePoints)
                    copyto!(unconf.TPointsₚ, unconf.TPoints)
                    copyto!(boundary.TPointsₚ, boundary.TPoints)
                end
            end

            # Recalculate interpolated proposed t-T path
            TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], model.ageSteps)

             # Calculate model ages for each grain
             anneal!(pr, Teq, model.dt, model.tSteps, TSteps, ZRDAAM())
             for i=1:length(zircons)
                 first_index = 1 + floor(Int64,(model.tInit - data.CrystAge[i])/model.dt)
                 @views CalcHeAgesₚ[i] = HeAgeSpherical(zircons[i],TSteps[first_index:end], pr[first_index:end,first_index:end], model)
             end

            # Calculate log likelihood of proposal
            AnnealedSigma .= simannealsigma.(n, data.HeAge_sigma; simannealmodel=model)
            llₚ = normpdf_ll(data.HeAge, AnnealedSigma, CalcHeAgesₚ)
            llₗ = normpdf_ll(data.HeAge, AnnealedSigma, CalcHeAges) # Recalulate last one too with new AnnealedSigma
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
                copyto!(CalcHeAges, CalcHeAgesₚ)

                # These are saved for ouput, but not critical to the function of the MCMC loop
                copyto!(TSteps, TStepsₚ)
                acceptanceDist[n] = true
            end

            # Record results for analysis and troubleshooting
            llDist[n] = normpdf_ll(data.HeAge, UnAnnealedSigma, CalcHeAges) # Recalculated to constant baseline
            nDist[n] = nPoints # Distribution of # of points
            HeAgeDist[:,n] .= CalcHeAges # Distribution of He ages

            # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
            TStepDist[:,n] .= TSteps # Distribution of T paths

        end
        return (TStepDist, HeAgeDist, nDist, llDist, acceptanceDist)
    end
    export MCMC_vartcryst

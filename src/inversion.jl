
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
    function simannealsigma(n::Integer, sigma_analytical::AbstractFloat; simannealmodel::NamedTuple=(σModel=25.0, σAnnealing=35.0, λAnnealing=10/10^5))
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

    # Check if point k is distinct from other points in list within ± δ
    function isdistinct(points::DenseArray, npoints::Int, k::Int, δ::Number)
        @inbounds for i = 1:npoints
            if i!=k && abs(points[i] - points[k]) < δ
                return false
            end
        end
        return true
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
    function MCMC_vartcryst(data::NamedTuple, model::NamedTuple, nPoints::Int, agePoints::DenseVector{T}, TPoints::DenseVector{T}, unconf::NamedTuple, boundary::NamedTuple) where T <: AbstractFloat
        # Sanitize inputs
        @assert firstindex(agePoints) === 1
        @assert firstindex(TPoints) === 1
        halfwidth = T.(data.halfwidth)::DenseVector{T}
        halfwidth = T.(data.halfwidth)::DenseVector{T}
        CrystAge = T.(data.CrystAge)::DenseVector{T}
        HeAge = T.(data.HeAge)::DenseVector{T}
        HeAge_sigma = T.(data.HeAge_sigma)::DenseVector{T}
        U = T.(data.U)::DenseVector{T}
        Th = T.(data.Th)::DenseVector{T}
        nsteps = model.nsteps::Int
        maxPoints = model.maxPoints::Int
        minPoints = (haskey(model, :minPoints) ? model.minPoints : 1)::Int
        simplified = model.simplified::Bool
        dTmax = T(model.dTmax)::T
        ageSteps = T.(model.ageSteps)::DenseVector{T}
        tSteps = T.(model.tSteps)::DenseVector{T}
        tInit = T(model.tInit)::T
        TInit = T(model.TInit)::T
        TNow = T(model.TNow)::T
        dt = T(model.dt)::T
        dr = T(model.dr)::T

        # Calculate number of boundary and unconformity points and allocate buffer for interpolating
        extraPoints = length(boundary.agePoints) + length(unconf.agePoints)
        agePointBuffer = similar(agePoints, maxPoints+extraPoints)::DenseVector{T}
        TPointBuffer = similar(agePoints, maxPoints+extraPoints)::DenseVector{T}
        knot_index = similar(ageSteps, Int)::DenseVector{Int}

        # Calculate model ages for initial proposal
        ages = collectto!(agePointBuffer, view(agePoints, 1:nPoints), boundary.agePoints, unconf.agePoints)::StridedVector{T}
        temperatures = collectto!(TPointBuffer, view(TPoints, 1:nPoints), boundary.TPoints, unconf.TPoints)::StridedVector{T}
        TSteps = linterp1s(ages, temperatures, ageSteps)::DenseVector{T}
        calcHeAges = similar(HeAge)::DenseVector{T}
        pr, Teq = anneal(dt, tSteps, TSteps, ZRDAAM()) # Damage annealing history
        pr::DenseMatrix{T}
        Teq::DenseVector{T}

        # Prepare a Mineral object for each analysis
        diffusionmodel = (DzEa=T(model.DzEa), DzD0=T(model.DzD0), DN17Ea=T(model.DN17Ea), DN17D0=T(model.DN17D0))
        zircons = Array{Zircon{T}}(undef, length(halfwidth))::Vector{Zircon{T}}
        for i=1:length(zircons)
            # Iterate through each grain, calculate the modeled age for each
            first_index = 1 + floor(Int64,(tInit - CrystAge[i])/dt)
            zircons[i] = Zircon(halfwidth[i], dr, U[i], Th[i], dt, ageSteps[first_index:end])
            calcHeAges[i] = HeAgeSpherical(zircons[i], @views(TSteps[first_index:end]), @views(pr[first_index:end,first_index:end]), diffusionmodel)::T
        end

        # Simulated annealing of uncertainty
        simannealmodel = (σModel=T(model.σModel), σAnnealing=T(model.σAnnealing), λAnnealing=T(model.λAnnealing))
        σₐ = simannealsigma.(1, HeAge_sigma; simannealmodel)::DenseVector{T}
        σ = sqrt.(T(model.σModel)^2 .+ HeAge_sigma.^2)

        # Log-likelihood for initial proposal
        ll = normpdf_ll(HeAge, σₐ, calcHeAges)
        if simplified
            ll -= log(nPoints)
        end

        # Variables to hold proposals
        llₚ = ll
        nPointsₚ = nPoints
        agePointsₚ = similar(agePoints)::DenseVector{T}
        TPointsₚ = similar(TPoints)::DenseVector{T}
        calcHeAgesₚ = similar(calcHeAges)::DenseVector{T}
        TStepsₚ = similar(TSteps)::DenseVector{T}

        # distributions to populate
        HeAgedist = Array{T}(undef, length(HeAge), nsteps)
        TStepdist = Array{T}(undef, length(tSteps), nsteps)
        lldist = Array{T}(undef, nsteps)
        ndist = zeros(Int, nsteps)
        acceptancedist = zeros(Bool, nsteps)

        # Standard deviations of Gaussian proposal ("jumping") distributions
        # for temperature and time
        σⱼt = tInit/60
        σⱼT = TInit/60
        k = 1 # Index of chosen t-T point

        # Proposal probabilities (must sum to 1)
        move = 0.64
        birth = 0.15 # Must equal death
        death = 0.15 # Must equal birth
        movebounds = 0.06

        # Number of times to attempt to satisfy reheating rate: # Should be large
        # enough that proposal probabilities are unchanged, but low enough to prevent
        # infinite loop
        nattempts = 100_000

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
                for attempt ∈ 1:nattempts

                # Move the age of one model point
                agePointsₚ[k] += randn() * σⱼt
                if agePointsₚ[k] < (0 + dt)
                    # Reflecting boundary condition at (0 + dt)
                    agePointsₚ[k] = (0 + dt) - (agePointsₚ[k] - (0 + dt))
                elseif agePointsₚ[k] > (tInit - dt)
                    # Reflecting boundary condition at (tInit - dt)
                    agePointsₚ[k] = (tInit - dt) - (agePointsₚ[k] - (tInit - dt))
                end

                # Move the Temperature of one model point
                TPointsₚ[k] += randn() * σⱼT
                if TPointsₚ[k] < TNow
                    # Reflecting boundary conditions at TNow (0)
                    TPointsₚ[k] = TNow - (TPointsₚ[k] - TNow)
                elseif TPointsₚ[k] > TInit
                    # Reflecting boundary conditions at TInit (0)
                    TPointsₚ[k] = TInit - (TPointsₚ[k] - TInit)
                end

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agePointBuffer, view(agePointsₚ, 1:nPointsₚ), boundary.agePoints, unconf.agePointsₚ)::StridedVector{T}
                temperatures = collectto!(TPointBuffer, view(TPointsₚ, 1:nPointsₚ), boundary.TPointsₚ, unconf.TPointsₚ)::StridedVector{T}
                linterp1s!(TStepsₚ, knot_index, ages, temperatures, ageSteps)

                # Retry unless we have satisfied the maximum reheating rate
                if isdistinct(agePointsₚ, nPointsₚ, k, dt) && maxdiff(TStepsₚ) < dTmax
                    break
                end
                # Copy last accepted solution to re-modify if we don't break
                copyto!(agePointsₚ, agePoints)
                copyto!(TPointsₚ, TPoints)
                if (attempt == nattempts)
                    @info """Warning: `move` proposals failed to satisfy reheating rate limit
                    maxdiff: $(maxdiff(TStepsₚ))
                    ages: $ages
                    temperatures: $temperatures"""
                end
                end
            elseif (r < move+birth) && (nPointsₚ < maxPoints)
                # Birth: add a new model point
                nPointsₚ += 1
                for attempt ∈ 1:nattempts
                agePointsₚ[nPointsₚ] = 0 + rand()*(tInit-0)
                TPointsₚ[nPointsₚ] = TNow + rand()*(TInit-TNow)

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agePointBuffer, view(agePointsₚ, 1:nPointsₚ), boundary.agePoints, unconf.agePointsₚ)::StridedVector{T}
                temperatures = collectto!(TPointBuffer, view(TPointsₚ, 1:nPointsₚ), boundary.TPointsₚ, unconf.TPointsₚ)::StridedVector{T}
                linterp1s!(TStepsₚ, knot_index, ages, temperatures, ageSteps)

                # Retry unless we have satisfied the maximum reheating rate
                if isdistinct(agePointsₚ, nPointsₚ, k, dt) && maxdiff(TStepsₚ) < dTmax
                    break
                end
                if (attempt == nattempts)
                    @info """Warning: new point proposals failed to satisfy reheating rate limit
                    maxdiff: $(maxdiff(TStepsₚ))
                    ages: $ages
                    temperatures: $temperatures
                    σⱼt: $(σⱼt)
                    σⱼT: $(σⱼT)"""
                end
                end
            elseif (r < move+birth+death) && (r >= move+birth) && (nPointsₚ > minPoints)
                # Death: remove a model point
                nPointsₚ -= 1 # Delete last point in array from proposal
                for attempt ∈ 1:nattempts
                k = ceil(Int, rand()*nPoints) # Choose point to delete
                agePointsₚ[k] = agePointsₚ[nPoints]
                TPointsₚ[k] = TPointsₚ[nPoints]

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agePointBuffer, view(agePointsₚ, 1:nPointsₚ), boundary.agePoints, unconf.agePointsₚ)::StridedVector{T}
                temperatures = collectto!(TPointBuffer, view(TPointsₚ, 1:nPointsₚ), boundary.TPointsₚ, unconf.TPointsₚ)::StridedVector{T}
                linterp1s!(TStepsₚ, knot_index, ages, temperatures, ageSteps)

                # Retry unless we have satisfied the maximum reheating rate
                (maxdiff(TStepsₚ) < dTmax) && break
                if (attempt == nattempts)
                    @info """Warning: point removal proposals failed to satisfy reheating rate limit
                    maxdiff: $(maxdiff(TStepsₚ))
                    ages: $ages
                    temperatures: $temperatures"""
                end
                end
            else
                # Move boundary conditions
                for attempt ∈ 1:nattempts
                # Move the temperatures of the starting or ending boundaries
                # @. boundary.TPointsₚ = boundary.T₀ + rand()*boundary.ΔT
                k = rand(1:length(boundary.TPointsₚ))
                boundary.TPointsₚ[k] = boundary.T₀[k] + rand()*boundary.ΔT[k]

                # If there's an imposed unconformity, adjust within parameters
                if length(unconf.agePoints) > 0
                    @. unconf.agePointsₚ = unconf.Age₀ + rand()*unconf.ΔAge
                    @. unconf.TPointsₚ = unconf.T₀ + rand()*unconf.ΔT
                end

                # Recalculate interpolated proposed t-T path
                ages = collectto!(agePointBuffer, view(agePointsₚ, 1:nPointsₚ), boundary.agePoints, unconf.agePointsₚ)::StridedVector{T}
                temperatures = collectto!(TPointBuffer, view(TPointsₚ, 1:nPointsₚ), boundary.TPointsₚ, unconf.TPointsₚ)::StridedVector{T}
                linterp1s!(TStepsₚ, knot_index, ages, temperatures, ageSteps)

                # Retry unless we have satisfied the maximum reheating rate
                (maxdiff(TStepsₚ) < dTmax) && break

                # Copy last accepted solution to re-modify if we don't break
                copyto!(unconf.agePointsₚ, unconf.agePoints)
                copyto!(unconf.TPointsₚ, unconf.TPoints)
                copyto!(boundary.TPointsₚ, boundary.TPoints)
                if (attempt == nattempts)
                    @info """Warning: `movebounds` proposals failed to satisfy reheating rate limit
                    maxdiff: $(maxdiff(TStepsₚ))
                    ages: $ages
                    temperatures: $temperatures"""
                end
                end
            end

            # ages = collectto!(agePointBuffer, view(agePointsₚ, 1:nPointsₚ), boundary.agePoints, unconf.agePointsₚ)::StridedVector{T}
            # temperatures = collectto!(TPointBuffer, view(TPointsₚ, 1:nPointsₚ), boundary.TPointsₚ, unconf.TPointsₚ)::StridedVector{T}
            # linterp1s!(TStepsₚ, knot_index, ages, temperatures, ageSteps)

             # Calculate model ages for each grain
            anneal!(pr, Teq, dt, tSteps, TStepsₚ, ZRDAAM())
            for i=1:length(zircons)
                first_index = 1 + floor(Int64,(tInit - CrystAge[i])/dt)
                calcHeAgesₚ[i] = HeAgeSpherical(zircons[i], @views(TStepsₚ[first_index:end]), @views(pr[first_index:end,first_index:end]), diffusionmodel)::T
            end

            # Calculate log likelihood of proposal
            σₐ .= simannealsigma.(n, HeAge_sigma; simannealmodel)
            llₚ = normpdf_ll(HeAge, σₐ, calcHeAgesₚ)
            llₗ = normpdf_ll(HeAge, σₐ, calcHeAges) # Recalulate last one too with new σₐ
            if simplified # slightly penalize more complex t-T paths
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

                # Update jumping distribution based on size of last accepted move
                if r < move
                    if agePointsₚ[k] != agePoints[k]
                        σⱼt = 3 * abs(agePointsₚ[k] - agePoints[k])
                    end
                    if TPointsₚ[k] != TPoints[k]
                        σⱼT = 3 * abs(TPointsₚ[k] - TPoints[k])
                    end
                end

                # Not critical to the function of the MCMC loop, but critical for recording stationary distribution!
                copyto!(TSteps, TStepsₚ)
                acceptancedist[n] = true
            end

            # Record results for analysis and troubleshooting
            lldist[n] = normpdf_ll(HeAge, σ, calcHeAges) # Recalculated to constant baseline
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

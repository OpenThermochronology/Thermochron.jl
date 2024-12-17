
    function modelages!(μcalc::AbstractVector{T}, σcalc::AbstractVector{T}, data::Vector{<:Chronometer{T}}, Tsteps::AbstractVector{T}, zdm::ZirconHeliumModel{T}, adm::ApatiteHeliumModel{T}, aftm::AnnealingModel{T}) where {T<:AbstractFloat}
        @assert eachindex(μcalc) == eachindex(σcalc) == eachindex(data)
        imax = argmax(i->length(data[i].agesteps), eachindex(data))
        tsteps = data[imax].tsteps
        tmax = last(tsteps)
        dt = step(tsteps)
        @assert issorted(tsteps)
        @assert eachindex(tsteps) == eachindex(Tsteps)

        # Pre-anneal ZRDAAM samples, if any
        isa(zdm, ZRDAAM) && anneal!(data, ZirconHe{T}, tsteps, Tsteps, zdm)
        # Pre-anneal RDAAM samples, if any
        isa(adm, RDAAM) && anneal!(data, ApatiteHe{T}, tsteps, Tsteps, adm)

        for i in eachindex(data)
            c = data[i]
            first_index = 1 + Int((tmax - last(c.tsteps))÷dt)
            if isa(c, ZirconHe)
                μcalc[i] = modelage(data[i], @views(Tsteps[first_index:end]), zdm)
            elseif isa(c, ApatiteHe)
                μcalc[i] = modelage(c, @views(Tsteps[first_index:end]), adm)
            elseif isa(c, ApatiteFT)
                μcalc[i] = modelage(c, @views(Tsteps[first_index:end]), aftm)
            elseif isa(c, ApatiteTrackLength)
                l,σ = modellength(c, @views(Tsteps[first_index:end]), aftm) .* aftm.l0
                μcalc[i] = l
                σcalc[i] = sqrt(σ^2 + aftm.l0_sigma^2)
            else
                # NaN if not calculated
                μcalc[i] = T(NaN)

            end
        end
        return μcalc, σcalc
    end

    """
    ```julia
    MCMC(data::Vector{<:Chronometer}, model::NamedTuple, npoints::Int, agepoints::Vector, Tpoints::Vector, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])
    ```
    Markov chain Monte Carlo time-Temperature inversion of the thermochronometric data 
    specified as a vector `chrons` of `Chronometer` objects (`ZirconHe`, `ApatiteHe`, 
    `ApatiteFT`, etc.) and model parameters specified by the named tuple `model`, 
    with variable diffusion kinetics.

    Returns a `TTResult` object containing posterior time-temperature paths,

    See also `MCMC_varkinetics` for a variant with variable diffusion kinetics.

    ## Examples
    ```julia
    tT = MCMC(data, model, npoints, agepoints, Tpoints, constraint, boundary)
    ```
    """
    MCMC(data::NamedTuple, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where {T <: AbstractFloat} = MCMC(chronometers(T, data, model), model, boundary, constraint, detail)
    function MCMC(data::Vector{<:ChronometerUnion{T}}, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where T <: AbstractFloat
        # Process inputs
        observed = val.(data)::Vector{T}
        observed_sigma = err.(data)::Vector{T}
        burnin = (haskey(model, :burnin) ? model.burnin : 5*10^5)::Int
        nsteps = (haskey(model, :nsteps) ? model.nsteps : 10^6)::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        maxpoints = (haskey(model, :maxpoints) ? model.maxpoints : 50)::Int
        npoints = (haskey(model, :npoints) ? model.npoints : minpoints)::Int
        totalpoints = maxpoints + boundary.npoints + constraint.npoints::Int
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicsigma = (haskey(model, :dynamicsigma) ? model.dynamicsigma : true)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : 5)::T
        agesteps = floatrange(model.agesteps)
        tsteps = floatrange(model.tsteps)
        @assert tsteps == reverse(agesteps)
        @assert issorted(tsteps)
        tnow, tinit = extrema(boundary.agepoints)
        Tnow, Tinit = extrema(boundary.T₀)
        Tr = T(haskey(model, :Tr) ? model.Tr : (Tinit+Tnow)/2)::T
        dt = T(model.dt)::T
        σmodel = T(haskey(model, :σmodel) ? model.σmodel : 0)::T
        σannealing = T(haskey(model, :σannealing) ? model.σannealing : 125)::T
        λannealing = T(haskey(model, :λannealing) ? model.λannealing : 2/burnin)::T

        # Arrays to hold all t and T points (up to npoints=maxpoints)
        agepoints = zeros(T, maxpoints) 
        Tpoints = zeros(T, maxpoints)
    
        # Fill some intermediate points to give the MCMC something to work with
        agepoints[1:npoints] .= range(tnow, tinit, length=npoints)
        Tpoints[1:npoints] .= Tr # Degrees C

        # Calculate number of boundary and unconformity points and allocate buffer for interpolating
        agepointbuffer = similar(agepoints, totalpoints)::Vector{T}
        Tpointbuffer = similar(agepoints, totalpoints)::Vector{T}
        knot_index = similar(agesteps, Int)::Vector{Int}

        # Prepare to calculate model ages for initial proposal
        ages = collectto!(agepointbuffer, view(agepoints, 1:npoints), boundary.agepoints, constraint.agepoints)::StridedVector{T}
        temperatures = collectto!(Tpointbuffer, view(Tpoints, 1:npoints), boundary.Tpoints, constraint.Tpoints)::StridedVector{T}
        Tsteps = linterp1s(ages, temperatures, agesteps)::Vector{T}
        μcalc = zeros(T, length(observed))
        σcalc = fill(T(σmodel), length(observed))

        # Damage models for each mineral
        zdm = (haskey(model, :zdm) ? model.zdm : ZRDAAM())::ZirconHeliumModel{T}
        adm = (haskey(model, :adm) ? model.adm : RDAAM())::ApatiteHeliumModel{T}
        aftm = (haskey(model, :aftm) ? model.aftm : FCKetcham2007)::AnnealingModel{T}

        # See what minerals we have
        (haszhe = any(x->isa(x, ZirconHe), data)) && @info "Inverting for He ages of $(count(x->isa(x, ZirconHe), data)) zircons"
        (hasahe = any(x->isa(x, ApatiteHe), data)) && @info "Inverting for He ages of $(count(x->isa(x, ApatiteHe), data)) apatites"
        (hasaft = any(x->isa(x, ApatiteFT), data)) && @info "Inverting for fission track ages of $(count(x->isa(x, ApatiteFT), data)) apatites"
        (hasatl = any(x->isa(x, ApatiteTrackLength), data)) && @info "Inverting for track lengths of $(count(x->isa(x, ApatiteTrackLength), data)) apatite fission tracks"
        
        # Standard deviations of Gaussian proposal ("jumping") distributions
        # for temperature and time
        σⱼt = fill((tinit-tnow)/60, maxpoints)
        σⱼT = fill((Tinit-Tnow)/60, maxpoints)

        # Simulated annealing of uncertainty
        σₐ = simannealsigma.(1, observed_sigma, σannealing, λannealing)::Vector{T}
        σ = observed_sigma

        # Log-likelihood for initial proposal
        modelages!(μcalc, σcalc, data, Tsteps, zdm, adm, aftm)
        llna = llnaₚ = diff_ll(Tsteps, dTmax, dTmax_sigma) + (simplified ? -log(npoints) : zero(T))
        ll = llₚ =  norm_ll(observed, σₐ, μcalc, σcalc) + llna

        # Variables to hold proposals
        npointsₚ = npoints
        agepointsₚ = copy(agepoints)::Vector{T}
        Tpointsₚ = copy(Tpoints)::Vector{T}
        μcalcₚ = copy(μcalc)::Vector{T}
        σcalcₚ = copy(σcalc)::Vector{T}
        σⱼtₚ = copy(σⱼt)
        σⱼTₚ = copy(σⱼT)

        # Proposal probabilities (must sum to 1)
        p_move = 0.64
        p_birth = 0.15 # Must equal p_death
        p_death = 0.15 # Must equal p_birth
        p_bounds = 0.06

        bprogress = Progress(burnin, dt=1, desc="MCMC burn-in ($(burnin) steps)")
        progress_interval = ceil(Int,sqrt(burnin))
        for n = 1:burnin
            if detail.minpoints > 0
                enoughpoints = min(pointsininterval(agepoints, npoints, detail.agemin, detail.agemax, dt), detail.minpoints)::Int
            end
            @label brestart

            # Copy proposal from last accepted solution
            npointsₚ = npoints
            copyto!(agepointsₚ, agepoints)
            copyto!(Tpointsₚ, Tpoints)
            copyto!(constraint.agepointsₚ, constraint.agepoints)
            copyto!(constraint.Tpointsₚ, constraint.Tpoints)
            copyto!(boundary.Tpointsₚ, boundary.Tpoints)
            copyto!(σcalcₚ, σcalc)
            copyto!(σⱼtₚ, σⱼt)
            copyto!(σⱼTₚ, σⱼT)

            # Randomly choose an option and point (if applicable) to adjust
            r = rand()
            k = rand(1:npoints)

            # Adjust the proposal
            if r < p_move
                # Move one t-T point
                movepoint!(agepointsₚ, Tpointsₚ, k, σⱼtₚ[k], σⱼTₚ[k], boundary)

            elseif (r < p_move+p_birth) && (npoints < maxpoints)
                # Birth: add a new model point
                k = npointsₚ = npoints + 1
                addpoint!(agepointsₚ, Tpointsₚ, σⱼtₚ, σⱼTₚ, k, boundary)

            elseif (r < p_move+p_birth+p_death) && (r >= p_move+p_birth) && (npoints > max(minpoints, detail.minpoints))
                # Death: remove a model point
                npointsₚ = npoints - 1
                agepointsₚ[k] = agepointsₚ[npoints]
                Tpointsₚ[k] = Tpointsₚ[npoints]
                σⱼtₚ[k] = σⱼtₚ[npoints]
                σⱼTₚ[k] = σⱼTₚ[npoints]

            elseif (r < p_move+p_birth+p_death+p_bounds)
                # Move the temperatures of the starting and ending boundaries
                movebounds!(boundary)
                # If there's an imposed unconformity or other t-T constraint, adjust within bounds
                movebounds!(constraint, boundary)

                # Move uncertainties
                dynamicsigma && movesigma!(σcalcₚ, data)

            end

            # Recalculate interpolated proposed t-T path
            ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, constraint.agepointsₚ)::StridedVector{T}
            temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, constraint.Tpointsₚ)::StridedVector{T}
            linterp1s!(Tsteps, knot_index, ages, temperatures, agesteps)

            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (pointsininterval(agepointsₚ, npointsₚ, detail.agemin, detail.agemax, dt) < enoughpoints) && @goto brestart
            end

            # Calculate model ages for each grain
            modelages!(μcalcₚ, σcalcₚ, data, Tsteps, zdm, adm, aftm)

            # Calculate log likelihood of proposal
            llnaₚ = diff_ll(Tsteps, dTmax, dTmax_sigma)
            simplified && (llnaₚ += -log(npointsₚ))
            σₐ .= simannealsigma.(n, observed_sigma, σannealing, λannealing)
            llₚ = norm_ll(observed, σₐ, μcalcₚ, σcalcₚ) + llnaₚ
            llₗ = norm_ll(observed, σₐ, μcalc, σcalc) + llna # Recalulate last one too with new σₐ

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - llₗ)

                # Update jumping distribution based on size of current accepted p_move
                if dynamicjumping && r < p_move
                    if agepointsₚ[k] != agepoints[k]
                        σⱼtₚ[k] = ℯ * abs(agepointsₚ[k] - agepoints[k])
                    end
                    if Tpointsₚ[k] != Tpoints[k]
                        σⱼTₚ[k] = ℯ * abs(Tpointsₚ[k] - Tpoints[k])
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
                copyto!(μcalc, μcalcₚ)
                copyto!(σcalc, σcalcₚ)
                copyto!(σⱼt, σⱼtₚ)
                copyto!(σⱼT, σⱼTₚ)
            end

            # Update progress meter every `progress_interval` steps
            (mod(n, progress_interval) == 0) && update!(bprogress, n)
        end
        finish!(bprogress)

        # Final log likelihood
        ll = norm_ll(observed, σ, μcalc, σcalc) + llna

        # distributions to populate
        tpointdist = fill(T(NaN), totalpoints, nsteps)
        Tpointdist = fill(T(NaN), totalpoints, nsteps)
        resultdist = fill(T(NaN), length(observed), nsteps)
        σⱼtdist = zeros(T, nsteps)
        σⱼTdist = zeros(T, nsteps)
        lldist = zeros(T, nsteps)
        ndist = zeros(Int, nsteps)
        acceptancedist = falses(nsteps)

        progress = Progress(nsteps, dt=1, desc="MCMC collection ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps
            if detail.minpoints > 0
                enoughpoints = min(pointsininterval(agepoints, npoints, detail.agemin, detail.agemax, dt), detail.minpoints)::Int
            end
            @label crestart

            # Copy proposal from last accepted solution
            npointsₚ = npoints
            copyto!(agepointsₚ, agepoints)
            copyto!(Tpointsₚ, Tpoints)
            copyto!(constraint.agepointsₚ, constraint.agepoints)
            copyto!(constraint.Tpointsₚ, constraint.Tpoints)
            copyto!(boundary.Tpointsₚ, boundary.Tpoints)
            copyto!(σcalcₚ, σcalc)
            copyto!(σⱼtₚ, σⱼt)
            copyto!(σⱼTₚ, σⱼT)

            # Randomly choose an option and point (if applicable) to adjust
            r = rand()
            k = rand(1:npoints)

            # Adjust the proposal
            if r < p_move
                # Move one t-T point
                movepoint!(agepointsₚ, Tpointsₚ, k, σⱼtₚ[k], σⱼTₚ[k], boundary)

            elseif (r < p_move+p_birth) && (npoints < maxpoints)
                # Birth: add a new model point
                k = npointsₚ = npoints + 1
                addpoint!(agepointsₚ, Tpointsₚ, σⱼtₚ, σⱼTₚ, k, boundary)

            elseif (r < p_move+p_birth+p_death) && (r >= p_move+p_birth) && (npoints > max(minpoints, detail.minpoints))
                # Death: remove a model point
                npointsₚ = npoints - 1
                agepointsₚ[k] = agepointsₚ[npoints]
                Tpointsₚ[k] = Tpointsₚ[npoints]
                σⱼtₚ[k] = σⱼtₚ[npoints]
                σⱼTₚ[k] = σⱼTₚ[npoints]

            elseif (r < p_move+p_birth+p_death+p_bounds)
                # Move the temperatures of the starting and ending boundaries
                movebounds!(boundary)
                # If there's an imposed unconformity or other t-T constraint, adjust within bounds
                movebounds!(constraint, boundary)

                # Move uncertainties
                dynamicsigma && movesigma!(σcalcₚ, data)

            end

            # Recalculate interpolated proposed t-T path
            ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, constraint.agepointsₚ)::StridedVector{T}
            temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, constraint.Tpointsₚ)::StridedVector{T}
            linterp1s!(Tsteps, knot_index, ages, temperatures, agesteps)
            
            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (pointsininterval(agepointsₚ, npointsₚ, detail.agemin, detail.agemax, dt) < enoughpoints) && @goto crestart
            end

            # Calculate model ages for each grain
            modelages!(μcalcₚ, σcalcₚ, data, Tsteps, zdm, adm, aftm)

            # Calculate log likelihood of proposal
            llnaₚ = diff_ll(Tsteps, dTmax, dTmax_sigma)
            simplified && (llnaₚ += -log(npointsₚ))
            llₚ = norm_ll(observed, σ, μcalcₚ, σcalcₚ) + llnaₚ

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - ll)

                # Update jumping distribution based on size of current accepted p_move
                if dynamicjumping && r < p_move
                    if agepointsₚ[k] != agepoints[k]
                        σⱼtₚ[k] = ℯ * abs(agepointsₚ[k] - agepoints[k])
                    end
                    if Tpointsₚ[k] != Tpoints[k]
                        σⱼTₚ[k] = ℯ * abs(Tpointsₚ[k] - Tpoints[k])
                    end
                end

                # Update the currently accepted proposal
                ll = llₚ
                npoints = npointsₚ
                copyto!(agepoints, 1, agepointsₚ, 1, npoints)
                copyto!(Tpoints, 1, Tpointsₚ, 1, npoints)
                copyto!(constraint.agepoints, constraint.agepointsₚ)
                copyto!(constraint.Tpoints, constraint.Tpointsₚ)
                copyto!(boundary.Tpoints, boundary.Tpointsₚ)
                copyto!(μcalc, μcalcₚ)
                copyto!(σcalc, σcalcₚ)
                copyto!(σⱼt, σⱼtₚ)
                copyto!(σⱼT, σⱼTₚ)

                # Not critical to the function of the MCMC loop, but critical for recording stationary distribution!
                acceptancedist[n] = true
            end

            # Record results for analysis and troubleshooting
            lldist[n] = ll
            ndist[n] = npoints # distribution of # of points
            σⱼtdist[n] = σⱼt[k]
            σⱼTdist[n] = σⱼT[k]
            resultdist[:,n] .= μcalc # distribution of He ages

            # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
            collectto!(view(tpointdist, :, n), view(agepoints, 1:npoints), boundary.agepoints, constraint.agepoints)
            collectto!(view(Tpointdist, :, n), view(Tpoints, 1:npoints), boundary.Tpoints, constraint.Tpoints)

            # Update progress meter every `progress_interval` steps
            (mod(n, progress_interval) == 0) && update!(progress, n)
        end
        finish!(progress)

        ttresult = TTResult(
            tpointdist, 
            Tpointdist, 
            ndist, 
            resultdist, 
            σⱼtdist, 
            σⱼTdist, 
            lldist, 
            acceptancedist,
        )
        return ttresult
    end


    """
    ```julia
    MCMC_varkinetics(chrons::Vector{<:Chronometer}, model::NamedTuple, npoints::Int, agepoints::Vector, Tpoints::Vector, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])
    ```
    Markov chain Monte Carlo time-Temperature inversion of the thermochronometric data 
    specified as a vector `chrons` of `Chronometer` objects (`ZirconHe`, `ApatiteHe`, 
    `ApatiteFT`, etc.) and model parameters specified by the named tuple `model`, 
    with variable diffusion kinetics.

    Returns a `TTResult` object containing posterior time-temperature paths,
    and a `KineticResult` object containing the posterior kinetic parameters.

    See also `MCMC` for a variant with constant diffusion kinetics.

    ## Examples
    ```julia
    tT, kinetics = MCMC_varkinetics(data, model, npoints, agepoints, Tpoints, constraint, boundary)
    ```
    """
    MCMC_varkinetics(data::NamedTuple, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where {T <: AbstractFloat} = MCMC_varkinetics(chronometers(T, data, model), model, boundary, constraint, detail)
    function MCMC_varkinetics(data::Vector{<:ChronometerUnion{T}}, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where T <: AbstractFloat
        # Process inputs
        observed = val.(data)::Vector{T}
        observed_sigma = err.(data)::Vector{T}
        burnin = (haskey(model, :burnin) ? model.burnin : 5*10^5)::Int
        nsteps = (haskey(model, :nsteps) ? model.nsteps : 10^6)::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        maxpoints = (haskey(model, :maxpoints) ? model.maxpoints : 50)::Int
        npoints = (haskey(model, :npoints) ? model.npoints : minpoints)::Int
        totalpoints = maxpoints + boundary.npoints + constraint.npoints::Int
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicsigma = (haskey(model, :dynamicsigma) ? model.dynamicsigma : true)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : 5)::T
        agesteps = floatrange(model.agesteps)
        tsteps = floatrange(model.tsteps)
        @assert tsteps == reverse(agesteps)
        @assert issorted(tsteps)
        tnow, tinit = extrema(boundary.agepoints)
        Tnow, Tinit = extrema(boundary.T₀)
        Tr = T(haskey(model, :Tr) ? model.Tr : (Tinit+Tnow)/2)::T
        dt = T(model.dt)::T
        σmodel = T(haskey(model, :σmodel) ? model.σmodel : 0)::T
        σannealing = T(haskey(model, :σannealing) ? model.σannealing : 125)::T
        λannealing = T(haskey(model, :λannealing) ? model.λannealing : 2/burnin)::T

        # Arrays to hold all t and T points (up to npoints=maxpoints)
        agepoints = zeros(T, maxpoints) 
        Tpoints = zeros(T, maxpoints)
    
        # Fill some intermediate points to give the MCMC something to work with
        agepoints[1:npoints] .= range(tnow, tinit, length=npoints)
        Tpoints[1:npoints] .= Tr # Degrees C

        # Calculate number of boundary and unconformity points and allocate buffer for interpolating
        agepointbuffer = similar(agepoints, totalpoints)::Vector{T}
        Tpointbuffer = similar(agepoints, totalpoints)::Vector{T}
        knot_index = similar(agesteps, Int)::Vector{Int}

        # Prepare to calculate model ages for initial proposal
        ages = collectto!(agepointbuffer, view(agepoints, 1:npoints), boundary.agepoints, constraint.agepoints)::StridedVector{T}
        temperatures = collectto!(Tpointbuffer, view(Tpoints, 1:npoints), boundary.Tpoints, constraint.Tpoints)::StridedVector{T}
        Tsteps = linterp1s(ages, temperatures, agesteps)::Vector{T}
        μcalc = zeros(T, length(observed))
        σcalc = fill(T(σmodel), length(observed))
        
        # Damage models for each mineral
        zdm₀ = zdm = zdmₚ = (haskey(model, :zdm) ? model.zdm : ZRDAAM())::ZirconHeliumModel{T}
        adm₀ = adm = admₚ =  (haskey(model, :adm) ? model.adm : RDAAM())::ApatiteHeliumModel{T}
        aftm = (haskey(model, :aftm) ? model.aftm : FCKetcham2007)::AnnealingModel{T}

        # See what minerals we have
        (haszhe = any(x->isa(x, ZirconHe), data)) && @info "Inverting for He ages of $(count(x->isa(x, ZirconHe), data)) zircons"
        (hasahe = any(x->isa(x, ApatiteHe), data)) && @info "Inverting for He ages of $(count(x->isa(x, ApatiteHe), data)) apatites"
        (hasaft = any(x->isa(x, ApatiteFT), data)) && @info "Inverting for fission track ages of $(count(x->isa(x, ApatiteFT), data)) apatites"
        (hasatl = any(x->isa(x, ApatiteTrackLength), data)) && @info "Inverting for track lengths of $(count(x->isa(x, ApatiteTrackLength), data)) apatite fission tracks"

        # Standard deviations of Gaussian proposal ("jumping") distributions
        # for temperature and time
        σⱼt = fill((tinit-tnow)/60, maxpoints)
        σⱼT = fill((Tinit-Tnow)/60, maxpoints)

        # Simulated annealing of uncertainty
        σₐ = simannealsigma.(1, observed_sigma, σannealing, λannealing)::Vector{T}
        σ = observed_sigma
        
        # Log-likelihood for initial proposal
        modelages!(μcalc, σcalc, data, Tsteps, zdm, adm, aftm)
        llna = llnaₚ = diff_ll(Tsteps, dTmax, dTmax_sigma) + loglikelihood(admₚ, adm₀) + loglikelihood(zdmₚ, zdm₀) + (simplified ? -log(npoints) : zero(T))
        ll = llₚ =  norm_ll(observed, σₐ, μcalc, σcalc) + llna

        # Variables to hold proposals
        npointsₚ = npoints
        agepointsₚ = copy(agepoints)::Vector{T}
        Tpointsₚ = copy(Tpoints)::Vector{T}
        μcalcₚ = copy(μcalc)::Vector{T}
        σcalcₚ = copy(σcalc)
        σⱼtₚ = copy(σⱼt)
        σⱼTₚ = copy(σⱼT)

        # Proposal probabilities (must sum to 1)
        p_move = 0.6
        p_birth = 0.15 # Must equal p_death
        p_death = 0.15 # Must equal p_birth
        p_bounds = 0.06
        p_kinetics = 0.04

        bprogress = Progress(burnin, dt=1, desc="MCMC burn-in ($(burnin) steps)")
        progress_interval = ceil(Int,sqrt(burnin))
        for n = 1:burnin
            if detail.minpoints > 0
                enoughpoints = min(pointsininterval(agepoints, npoints, detail.agemin, detail.agemax, dt), detail.minpoints)::Int
            end
            @label brestart

            # Copy proposal from last accepted solution
            admₚ = adm
            zdmₚ = zdm
            npointsₚ = npoints
            copyto!(agepointsₚ, agepoints)
            copyto!(Tpointsₚ, Tpoints)
            copyto!(constraint.agepointsₚ, constraint.agepoints)
            copyto!(constraint.Tpointsₚ, constraint.Tpoints)
            copyto!(boundary.Tpointsₚ, boundary.Tpoints)
            copyto!(σcalcₚ, σcalc)
            copyto!(σⱼtₚ, σⱼt)
            copyto!(σⱼTₚ, σⱼT)

            # Randomly choose an option and point (if applicable) to adjust
            r = rand()
            k = rand(1:npoints)

            # Adjust the proposal
            if r < p_move
                # Move one t-T point
                movepoint!(agepointsₚ, Tpointsₚ, k, σⱼtₚ[k], σⱼTₚ[k], boundary)

            elseif (r < p_move+p_birth) && (npoints < maxpoints)
                # Birth: add a new model point
                k = npointsₚ = npoints + 1
                addpoint!(agepointsₚ, Tpointsₚ, σⱼtₚ, σⱼTₚ, k, boundary)

            elseif (r < p_move+p_birth+p_death) && (r >= p_move+p_birth) && (npoints > max(minpoints, detail.minpoints))
                # Death: remove a model point
                npointsₚ = npoints - 1
                agepointsₚ[k] = agepointsₚ[npoints]
                Tpointsₚ[k] = Tpointsₚ[npoints]
                σⱼtₚ[k] = σⱼtₚ[npoints]
                σⱼTₚ[k] = σⱼTₚ[npoints]

            elseif (r < p_move+p_birth+p_death+p_bounds)
                # Move the temperatures of the starting and ending boundaries
                movebounds!(boundary)
                # If there's an imposed unconformity or other t-T constraint, adjust within bounds
                movebounds!(constraint, boundary)

                # Move uncertainties
                dynamicsigma && movesigma!(σcalcₚ, data)

            elseif (r < p_move+p_birth+p_death+p_bounds+p_kinetics)
                # Adjust kinetic parameters, one at a time
                haszhe && (zdmₚ = movekinetics(zdm))
                hasahe && (admₚ = movekinetics(adm))

            end

            # Recalculate interpolated proposed t-T path
            ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, constraint.agepointsₚ)::StridedVector{T}
            temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, constraint.Tpointsₚ)::StridedVector{T}
            linterp1s!(Tsteps, knot_index, ages, temperatures, agesteps)

            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (pointsininterval(agepointsₚ, npointsₚ, detail.agemin, detail.agemax, dt) < enoughpoints) && @goto brestart
            end
               
            # Calculate model ages for each grain
            modelages!(μcalcₚ, σcalcₚ, data, Tsteps, zdmₚ, admₚ, aftm)

            # Calculate log likelihood of proposal
            llnaₚ = diff_ll(Tsteps, dTmax, dTmax_sigma)
            llnaₚ += loglikelihood(admₚ, adm₀) + loglikelihood(zdmₚ, zdm₀)
            simplified && (llnaₚ += -log(npointsₚ))
            σₐ .= simannealsigma.(n, observed_sigma, σannealing, λannealing)
            llₚ = norm_ll(observed, σₐ, μcalcₚ, σcalcₚ) + llnaₚ
            llₗ = norm_ll(observed, σₐ, μcalc, σcalc) + llna # Recalulate last one too with new σₐ

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - llₗ)

                # Update jumping distribution based on size of current accepted p_move
                if dynamicjumping && r < p_move
                    if agepointsₚ[k] != agepoints[k]
                        σⱼtₚ[k] = ℯ * abs(agepointsₚ[k] - agepoints[k])
                    end
                    if Tpointsₚ[k] != Tpoints[k]
                        σⱼTₚ[k] = ℯ * abs(Tpointsₚ[k] - Tpoints[k])
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
                copyto!(μcalc, μcalcₚ)
                copyto!(σcalc, σcalcₚ)
                copyto!(σⱼt, σⱼtₚ)
                copyto!(σⱼT, σⱼTₚ)
            end

            # Update progress meter every `progress_interval` steps
            (mod(n, progress_interval) == 0) && update!(bprogress, n)
        end
        finish!(bprogress)

        # Final log likelihood
        ll = norm_ll(observed, σ, μcalc, σcalc) + llna

        # distributions to populate
        tpointdist = fill(T(NaN), totalpoints, nsteps)
        Tpointdist = fill(T(NaN), totalpoints, nsteps)
        resultdist = fill(T(NaN), length(observed), nsteps)
        σⱼtdist = zeros(T, nsteps)
        σⱼTdist = zeros(T, nsteps)
        lldist = zeros(T, nsteps)
        ndist = zeros(Int, nsteps)
        acceptancedist = falses(nsteps)
        admdist = Array{typeof(adm)}(undef, nsteps)
        zdmdist = Array{typeof(zdm)}(undef, nsteps)

        progress = Progress(nsteps, dt=1, desc="MCMC collection ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps
            if detail.minpoints > 0
                enoughpoints = min(pointsininterval(agepoints, npoints, detail.agemin, detail.agemax, dt), detail.minpoints)::Int
            end
            @label crestart

            # Copy proposal from last accepted solution
            admₚ = adm
            zdmₚ = zdm
            npointsₚ = npoints
            copyto!(agepointsₚ, agepoints)
            copyto!(Tpointsₚ, Tpoints)
            copyto!(constraint.agepointsₚ, constraint.agepoints)
            copyto!(constraint.Tpointsₚ, constraint.Tpoints)
            copyto!(boundary.Tpointsₚ, boundary.Tpoints)
            copyto!(σcalcₚ, σcalc)
            copyto!(σⱼtₚ, σⱼt)
            copyto!(σⱼTₚ, σⱼT)

            # Randomly choose an option and point (if applicable) to adjust
            r = rand()
            k = rand(1:npoints)

            # Adjust the proposal
            if r < p_move
                # Move one t-T point
                movepoint!(agepointsₚ, Tpointsₚ, k, σⱼtₚ[k], σⱼTₚ[k], boundary)

            elseif (r < p_move+p_birth) && (npoints < maxpoints)
                # Birth: add a new model point
                k = npointsₚ = npoints + 1
                addpoint!(agepointsₚ, Tpointsₚ, σⱼtₚ, σⱼTₚ, k, boundary)

            elseif (r < p_move+p_birth+p_death) && (r >= p_move+p_birth) && (npoints > max(minpoints, detail.minpoints))
                # Death: remove a model point
                npointsₚ = npoints - 1
                agepointsₚ[k] = agepointsₚ[npoints]
                Tpointsₚ[k] = Tpointsₚ[npoints]
                σⱼtₚ[k] = σⱼtₚ[npoints]
                σⱼTₚ[k] = σⱼTₚ[npoints]

            elseif (r < p_move+p_birth+p_death+p_bounds)
                # Move the temperatures of the starting and ending boundaries
                movebounds!(boundary)
                # If there's an imposed unconformity or other t-T constraint, adjust within bounds
                movebounds!(constraint, boundary)

                # Move uncertainties
                dynamicsigma && movesigma!(σcalcₚ, data)

            elseif (r < p_move+p_birth+p_death+p_bounds+p_kinetics)
                # Adjust kinetic parameters, one at a time
                haszhe && (zdmₚ = movekinetics(zdm))
                hasahe && (admₚ = movekinetics(adm))

            end

            # Recalculate interpolated proposed t-T path
            ages = collectto!(agepointbuffer, view(agepointsₚ, 1:npointsₚ), boundary.agepoints, constraint.agepointsₚ)::StridedVector{T}
            temperatures = collectto!(Tpointbuffer, view(Tpointsₚ, 1:npointsₚ), boundary.Tpointsₚ, constraint.Tpointsₚ)::StridedVector{T}
            linterp1s!(Tsteps, knot_index, ages, temperatures, agesteps)

            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (pointsininterval(agepointsₚ, npointsₚ, detail.agemin, detail.agemax, dt) < enoughpoints) && @goto crestart
            end

            # Calculate model ages for each grain
            modelages!(μcalcₚ, σcalcₚ, data, Tsteps, zdmₚ, admₚ, aftm)

            # Calculate log likelihood of proposal
            llnaₚ = diff_ll(Tsteps, dTmax, dTmax_sigma)
            llnaₚ += loglikelihood(admₚ, adm₀) + loglikelihood(zdmₚ, zdm₀)
            simplified && (llnaₚ += -log(npointsₚ))
            llₚ = norm_ll(observed, σ, μcalcₚ, σcalcₚ) + llnaₚ

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - ll)

                # Update jumping distribution based on size of current accepted p_move
                if dynamicjumping && r < p_move
                    if agepointsₚ[k] != agepoints[k]
                        σⱼtₚ[k] = ℯ * abs(agepointsₚ[k] - agepoints[k])
                    end
                    if Tpointsₚ[k] != Tpoints[k]
                        σⱼTₚ[k] = ℯ * abs(Tpointsₚ[k] - Tpoints[k])
                    end
                end

                # Update the currently accepted proposal
                adm = admₚ
                zdm = zdmₚ
                ll = llₚ
                npoints = npointsₚ
                copyto!(agepoints, 1, agepointsₚ, 1, npoints)
                copyto!(Tpoints, 1, Tpointsₚ, 1, npoints)
                copyto!(constraint.agepoints, constraint.agepointsₚ)
                copyto!(constraint.Tpoints, constraint.Tpointsₚ)
                copyto!(boundary.Tpoints, boundary.Tpointsₚ)
                copyto!(μcalc, μcalcₚ)
                copyto!(σcalc, σcalcₚ)
                copyto!(σⱼt, σⱼtₚ)
                copyto!(σⱼT, σⱼTₚ)

                # Not critical to the function of the MCMC loop, but critical for recording stationary distribution!
                acceptancedist[n] = true
            end

            # Record results for analysis and troubleshooting
            lldist[n] = ll
            ndist[n] = npoints # distribution of # of points
            σⱼtdist[n] = σⱼt[k]
            σⱼTdist[n] = σⱼT[k]
            resultdist[:,n] .= μcalc # distribution of He ages
            admdist[n] = adm
            zdmdist[n] = zdm

            # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
            collectto!(view(tpointdist, :, n), view(agepoints, 1:npoints), boundary.agepoints, constraint.agepoints)
            collectto!(view(Tpointdist, :, n), view(Tpoints, 1:npoints), boundary.Tpoints, constraint.Tpoints)

            # Update progress meter every `progress_interval` steps
            (mod(n, progress_interval) == 0) && update!(progress, n)
        end
        finish!(progress)

        ttresult = TTResult(
            tpointdist, 
            Tpointdist, 
            ndist, 
            resultdist, 
            σⱼtdist, 
            σⱼTdist, 
            lldist, 
            acceptancedist,
        )
        kineticresult = KineticResult(
            admdist,
            zdmdist
        )
        return ttresult, kineticresult
    end

    """
    ```julia
    MCMC(data::Vector{<:Chronometer}, model::NamedTuple, npoints::Int, path.agepoints::Vector, path.Tpoints::Vector, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])
    ```
    Markov chain Monte Carlo time-Temperature inversion of the thermochronometric data 
    specified as a vector `chrons` of `Chronometer` objects (`ZirconHe`, `ApatiteHe`, 
    `ApatiteFT`, etc.) and model parameters specified by the named tuple `model`, 
    with variable diffusion kinetics.

    Returns a `TTResult` object containing posterior time-temperature paths,

    See also `MCMC_varkinetics` for a variant with variable diffusion kinetics.

    ## Examples
    ```julia
    tT = MCMC(data::NamedTuple, model::NamedTuple, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])
    ```
    """
    MCMC(data::NamedTuple, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where {T <: AbstractFloat} = MCMC(chronometers(T, data, model), model, boundary, constraint, detail)
    function MCMC(data::Vector{<:ChronometerUnion{T}}, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where T <: AbstractFloat
        # Process inputs
        burnin = (haskey(model, :burnin) ? model.burnin : 5*10^5)::Int
        nsteps = (haskey(model, :nsteps) ? model.nsteps : 10^6)::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        maxpoints = (haskey(model, :maxpoints) ? model.maxpoints : 50)::Int
        npoints = (haskey(model, :npoints) ? model.npoints : minpoints)::Int
        npoints = max(npoints, detail.minpoints+1)
        totalpoints = maxpoints + boundary.npoints + constraint.npoints::Int
        trackhist = (haskey(model, :trackhist) ? model.trackhist : false)::Bool
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicsigma = (haskey(model, :dynamicsigma) ? model.dynamicsigma : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : 5)::T
        agesteps = floatrange(model.agesteps)
        tsteps = floatrange(model.tsteps)
        @assert tsteps == reverse(agesteps)
        @assert issorted(tsteps)
        dt = T(model.dt)::T
        σmodel = T(haskey(model, :σmodel) ? model.σmodel : 0)::T
        T0annealing = T(haskey(model, :T0annealing) ? model.T0annealing : 1)::T
        λannealing = T(haskey(model, :λannealing) ? model.λannealing : 7/burnin)::T

        # Damage models for each mineral
        zdm = (haskey(model, :zdm) ? model.zdm : ZRDAAM())::ZirconHeliumModel{T}
        adm = (haskey(model, :adm) ? model.adm : RDAAM())::ApatiteHeliumModel{T}
        zftm = (haskey(model, :zftm) ? model.zftm : Yamada2007PC())::ZirconAnnealingModel{T}
        mftm = (haskey(model, :mftm) ? model.mftm : Jones2021FA())::MonaziteAnnealingModel{T}
        aftm = (haskey(model, :aftm) ? model.aftm : Ketcham2007FC())::ApatiteAnnealingModel{T}

        # See what minerals we have
        (haszhe = any(x->isa(x, ZirconHe), data)) && @info "Inverting for He ages of $(count(x->isa(x, ZirconHe), data)) zircons"
        (hasahe = any(x->isa(x, ApatiteHe), data)) && @info "Inverting for He ages of $(count(x->isa(x, ApatiteHe), data)) apatites"
        (any(x->isa(x, GenericHe), data)) && @info "Inverting for He ages of $(count(x->isa(x, GenericHe), data)) generic He chronometers"
        (any(x->isa(x, ZirconFT), data)) && @info "Inverting for fission track ages of $(count(x->isa(x, ZirconFT), data)) zircons"
        (any(x->isa(x, ApatiteFT), data)) && @info "Inverting for fission track ages of $(count(x->isa(x, ApatiteFT), data)) apatites"
        (any(x->isa(x, ApatiteTrackLength), data)) && @info "Inverting for track lengths of $(count(x->isa(x, ApatiteTrackLength), data)) apatite fission tracks"
        (any(x->isa(x, GenericAr), data)) && @info "Inverting for Ar ages of $(count(x->isa(x, GenericAr), data)) generic Ar chronometers"

        # Struct to hold t-T path proposals and related variables
        path = TtPath(agesteps, constraint, boundary, detail, maxpoints)

        # Initial propopsal
        initialproposal!(path, npoints, dTmax) 

        # Prepare to calculate model ages for initial proposal
        μcalc = zeros(T, length(data))
        σcalc = fill(T(σmodel), length(data))

        # Log-likelihood for initial proposal
        ll = llₚ = model!(μcalc, σcalc, data, path.Tsteps, zdm, adm, zftm, mftm, aftm; trackhist) + diff_ll(path.Tsteps, dTmax, dTmax_sigma) + 
            (simplified ? -log(npoints) : zero(T))  + (dynamicsigma ? sum(x->-log1p(x), σcalc) : zero(T)) 

        # Variables to hold proposals
        npointsₚ = npoints
        μcalcₚ = copy(μcalc)::Vector{T}
        σcalcₚ = copy(σcalc)::Vector{T}

        # Proposal probabilities (must sum to 1)
        p_move = 0.64
        p_birth = 0.15 # Must equal p_death
        p_death = 0.15 # Must equal p_birth
        p_bounds = 0.06

        bprogress = Progress(burnin, dt=1, desc="MCMC burn-in ($(burnin) steps)")
        progress_interval = ceil(Int,sqrt(burnin))
        for n = 1:burnin
            if detail.minpoints > 0
                enoughpoints = min(pointsininterval(path.agepoints, npoints, detail.agemin, detail.agemax, dt), detail.minpoints)::Int
            end
            @label brestart

            # Copy proposal from last accepted solution
            resetproposal!(path)
            npointsₚ = npoints
            copyto!(σcalcₚ, σcalc)

            # Randomly choose an option and point (if applicable) to adjust
            r = rand()
            k = rand(Base.OneTo(npoints))

            # Adjust the proposal
            if r < p_move
                # Move one t-T point
                movepoint!(path, k)

            elseif (r < p_move+p_birth) && (npoints < maxpoints)
                # Birth: add a new model point
                k = npointsₚ = npoints + 1
                addpoint!(path, k)

            elseif (r < p_move+p_birth+p_death) && (r >= p_move+p_birth) && (npoints > max(minpoints, detail.minpoints))
                # Death: remove a model point
                npointsₚ = npoints - 1
                replacepoint!(path, k, npoints)

            elseif (r < p_move+p_birth+p_death+p_bounds)
                # Move the temperatures of the starting and ending boundaries
                # If there's an imposed unconformity or other t-T constraint, adjust within bounds
                movebounds!(path)

                # Move uncertainties
                dynamicsigma && movesigma!(σcalcₚ, data)

            end

            # Recalculate interpolated proposed t-T path
            collectproposal!(path, npointsₚ)

            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (pointsininterval(path.agepointsₚ, npointsₚ, detail.agemin, detail.agemax, dt) < enoughpoints) && @goto brestart
            end

            # Calculate model ages for each grain, log likelihood of proposal
            llₚ = model!(μcalcₚ, σcalcₚ, data, path.Tsteps, zdm, adm, zftm, mftm, aftm; trackhist)
            llₚ += diff_ll(path.Tsteps, dTmax, dTmax_sigma)
            simplified && (llₚ += -log(npointsₚ))
            dynamicsigma && (llₚ += sum(x->-log1p(x), σcalcₚ)) 

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - ll) / simannealT(n, T0annealing, λannealing)

                # Update jumping distribution based on size of current accepted p_move
                if dynamicjumping && r < p_move
                    if path.agepointsₚ[k] != path.agepoints[k]
                        path.σⱼtₚ[k] = ℯ * abs(path.agepointsₚ[k] - path.agepoints[k])
                    end
                    if path.Tpointsₚ[k] != path.Tpoints[k]
                        path.σⱼTₚ[k] = ℯ * abs(path.Tpointsₚ[k] - path.Tpoints[k])
                    end
                end

                # Update the currently accepted proposal
                acceptproposal!(path)
                ll = llₚ
                npoints = npointsₚ
                copyto!(μcalc, μcalcₚ)
                copyto!(σcalc, σcalcₚ)
            end

            # Update progress meter every `progress_interval` steps
            (mod(n, progress_interval) == 0) && update!(bprogress, n)
        end
        finish!(bprogress)

        # distributions to populate
        tpointdist = fill(T(NaN), totalpoints, nsteps)
        Tpointdist = fill(T(NaN), totalpoints, nsteps)
        resultdist = fill(T(NaN), length(data), nsteps)
        σⱼtdist = zeros(T, nsteps)
        σⱼTdist = zeros(T, nsteps)
        lldist = zeros(T, nsteps)
        ndist = zeros(Int, nsteps)
        acceptancedist = falses(nsteps)

        progress = Progress(nsteps, dt=1, desc="MCMC collection ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps
            if detail.minpoints > 0
                enoughpoints = min(pointsininterval(path.agepoints, npoints, detail.agemin, detail.agemax, dt), detail.minpoints)::Int
            end
            @label crestart

            # Copy proposal from last accepted solution
            resetproposal!(path)
            npointsₚ = npoints
            copyto!(σcalcₚ, σcalc)

            # Randomly choose an option and point (if applicable) to adjust
            r = rand()
            k = rand(Base.OneTo(npoints))

            # Adjust the proposal
            if r < p_move
                # Move one t-T point
                movepoint!(path, k)

            elseif (r < p_move+p_birth) && (npoints < maxpoints)
                # Birth: add a new model point
                k = npointsₚ = npoints + 1
                addpoint!(path, k)

            elseif (r < p_move+p_birth+p_death) && (r >= p_move+p_birth) && (npoints > max(minpoints, detail.minpoints))
                # Death: remove a model point
                npointsₚ = npoints - 1
                replacepoint!(path, k, npoints)

            elseif (r < p_move+p_birth+p_death+p_bounds)
                # Move the temperatures of the starting and ending boundaries
                # If there's an imposed unconformity or other t-T constraint, adjust within bounds
                movebounds!(path)

                # Move uncertainties
                dynamicsigma && movesigma!(σcalcₚ, data)

            end

            # Recalculate interpolated proposed t-T path
            collectproposal!(path, npointsₚ)
            
            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (pointsininterval(path.agepointsₚ, npointsₚ, detail.agemin, detail.agemax, dt) < enoughpoints) && @goto crestart
            end

            # Calculate model ages for each grain, log likelihood of proposal
            llₚ = model!(μcalcₚ, σcalcₚ, data, path.Tsteps, zdm, adm, zftm, mftm, aftm; trackhist)
            llₚ += diff_ll(path.Tsteps, dTmax, dTmax_sigma)
            simplified && (llₚ += -log(npointsₚ))
            dynamicsigma && (llₚ += sum(x->-log1p(x), σcalcₚ)) 

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - ll)

                # Update jumping distribution based on size of current accepted p_move
                if dynamicjumping && r < p_move
                    if path.agepointsₚ[k] != path.agepoints[k]
                        path.σⱼtₚ[k] = ℯ * abs(path.agepointsₚ[k] - path.agepoints[k])
                    end
                    if path.Tpointsₚ[k] != path.Tpoints[k]
                        path.σⱼTₚ[k] = ℯ * abs(path.Tpointsₚ[k] - path.Tpoints[k])
                    end
                end

                # Update the currently accepted proposal
                acceptproposal!(path)
                ll = llₚ
                npoints = npointsₚ                
                copyto!(μcalc, μcalcₚ)
                copyto!(σcalc, σcalcₚ)

                # Not critical to the function of the MCMC loop, but critical for recording stationary distribution!
                acceptancedist[n] = true
            end

            # Record results for analysis and troubleshooting
            lldist[n] = ll
            ndist[n] = npoints # distribution of # of points
            σⱼtdist[n] = path.σⱼt[k]
            σⱼTdist[n] = path.σⱼT[k]
            resultdist[:,n] .= μcalc # distribution of He ages

            # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
            collectto!(view(tpointdist, :, n), view(path.agepoints, Base.OneTo(npoints)), path.boundary.agepoints, path.constraint.agepoints)
            collectto!(view(Tpointdist, :, n), view(path.Tpoints, Base.OneTo(npoints)), path.boundary.Tpoints, path.constraint.Tpoints)

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
    MCMC_varkinetics(chrons::Vector{<:Chronometer}, model::NamedTuple, npoints::Int, path.agepoints::Vector, path.Tpoints::Vector, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])
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
    tT, kinetics = MCMC_varkinetics(data::NamedTuple, model::NamedTuple, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])
    ```
    """
    MCMC_varkinetics(data::NamedTuple, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where {T <: AbstractFloat} = MCMC_varkinetics(chronometers(T, data, model), model, boundary, constraint, detail)
    function MCMC_varkinetics(data::Vector{<:ChronometerUnion{T}}, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where T <: AbstractFloat
        # Process inputs
        burnin = (haskey(model, :burnin) ? model.burnin : 5*10^5)::Int
        nsteps = (haskey(model, :nsteps) ? model.nsteps : 10^6)::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        maxpoints = (haskey(model, :maxpoints) ? model.maxpoints : 50)::Int
        npoints = (haskey(model, :npoints) ? model.npoints : minpoints)::Int
        npoints = max(npoints, detail.minpoints+1)
        totalpoints = maxpoints + boundary.npoints + constraint.npoints::Int
        trackhist = (haskey(model, :trackhist) ? model.trackhist : false)::Bool
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicsigma = (haskey(model, :dynamicsigma) ? model.dynamicsigma : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : 5)::T
        agesteps = floatrange(model.agesteps)
        tsteps = floatrange(model.tsteps)
        @assert tsteps == reverse(agesteps)
        @assert issorted(tsteps)
        dt = T(model.dt)::T
        σmodel = T(haskey(model, :σmodel) ? model.σmodel : 0)::T
        T0annealing = T(haskey(model, :T0annealing) ? model.T0annealing : 1)::T
        λannealing = T(haskey(model, :λannealing) ? model.λannealing : 7/burnin)::T
        
        # Damage models for each mineral
        zdm₀ = zdm = zdmₚ = (haskey(model, :zdm) ? model.zdm : ZRDAAM())::ZirconHeliumModel{T}
        adm₀ = adm = admₚ =  (haskey(model, :adm) ? model.adm : RDAAM())::ApatiteHeliumModel{T}
        zftm = (haskey(model, :zftm) ? model.zftm : Yamada2007PC())::ZirconAnnealingModel{T}
        mftm = (haskey(model, :mftm) ? model.mftm : Jones2021FA())::MonaziteAnnealingModel{T}
        aftm = (haskey(model, :aftm) ? model.aftm : Ketcham2007FC())::ApatiteAnnealingModel{T}

        # See what minerals we have
        (haszhe = any(x->isa(x, ZirconHe), data)) && @info "Inverting for He ages of $(count(x->isa(x, ZirconHe), data)) zircons"
        (hasahe = any(x->isa(x, ApatiteHe), data)) && @info "Inverting for He ages of $(count(x->isa(x, ApatiteHe), data)) apatites"
        (any(x->isa(x, GenericHe), data)) && @info "Inverting for He ages of $(count(x->isa(x, GenericHe), data)) generic He chronometers"
        (any(x->isa(x, ZirconFT), data)) && @info "Inverting for fission track ages of $(count(x->isa(x, ZirconFT), data)) zircons"
        (any(x->isa(x, ApatiteFT), data)) && @info "Inverting for fission track ages of $(count(x->isa(x, ApatiteFT), data)) apatites"
        (any(x->isa(x, ApatiteTrackLength), data)) && @info "Inverting for track lengths of $(count(x->isa(x, ApatiteTrackLength), data)) apatite fission tracks"
        (any(x->isa(x, GenericAr), data)) && @info "Inverting for Ar ages of $(count(x->isa(x, GenericAr), data)) generic Ar chronometers"
        
        # Struct to hold t-T path proposals and related variables
        path = TtPath(agesteps, constraint, boundary, detail, maxpoints)

        # Initial propopsal
        initialproposal!(path, npoints, dTmax)
        
        # Prepare to calculate model ages for initial proposal
        μcalc = zeros(T, length(data))
        σcalc = fill(T(σmodel), length(data))

        # Log-likelihood for initial proposal
        ll = llₚ = model!(μcalc, σcalc, data, path.Tsteps, zdm, adm, zftm, mftm, aftm; trackhist) + 
            diff_ll(path.Tsteps, dTmax, dTmax_sigma) + kinetic_ll(admₚ, adm₀) + kinetic_ll(zdmₚ, zdm₀) + 
            (simplified ? -log(npoints) : zero(T)) + (dynamicsigma ? sum(x->-log1p(x), σcalc) : zero(T)) 
        
        # Variables to hold proposals
        npointsₚ = npoints
        μcalcₚ = copy(μcalc)::Vector{T}
        σcalcₚ = copy(σcalc)::Vector{T}

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
                enoughpoints = min(pointsininterval(path.agepoints, npoints, detail.agemin, detail.agemax, dt), detail.minpoints)::Int
            end
            @label brestart

            # Copy proposal from last accepted solution
            resetproposal!(path)
            admₚ = adm
            zdmₚ = zdm
            npointsₚ = npoints
            copyto!(σcalcₚ, σcalc)

            # Randomly choose an option and point (if applicable) to adjust
            r = rand()
            k = rand(Base.OneTo(npoints))

            # Adjust the proposal
            if r < p_move
                # Move one t-T point
                movepoint!(path, k)

            elseif (r < p_move+p_birth) && (npoints < maxpoints)
                # Birth: add a new model point
                k = npointsₚ = npoints + 1
                addpoint!(path, k)

            elseif (r < p_move+p_birth+p_death) && (r >= p_move+p_birth) && (npoints > max(minpoints, detail.minpoints))
                # Death: remove a model point
                npointsₚ = npoints - 1
                replacepoint!(path, k, npoints)

            elseif (r < p_move+p_birth+p_death+p_bounds)
                # Move the temperatures of the starting and ending boundaries
                # If there's an imposed unconformity or other t-T constraint, adjust within bounds
                movebounds!(path)

                # Move uncertainties
                dynamicsigma && movesigma!(σcalcₚ, data)

            elseif (r < p_move+p_birth+p_death+p_bounds+p_kinetics)
                # Adjust kinetic parameters, one at a time
                haszhe && (zdmₚ = movekinetics(zdm))
                hasahe && (admₚ = movekinetics(adm))

            end

            # Recalculate interpolated proposed t-T path
            collectproposal!(path, npointsₚ)

            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (pointsininterval(path.agepointsₚ, npointsₚ, detail.agemin, detail.agemax, dt) < enoughpoints) && @goto brestart
            end
               
            # Calculate model ages for each grain, log likelihood of proposal
            llₚ = model!(μcalcₚ, σcalcₚ, data, path.Tsteps, zdmₚ, admₚ, zftm, mftm, aftm; trackhist)
            llₚ += diff_ll(path.Tsteps, dTmax, dTmax_sigma)
            llₚ += kinetic_ll(admₚ, adm₀) + kinetic_ll(zdmₚ, zdm₀)
            simplified && (llₚ += -log(npointsₚ))
            dynamicsigma && (llₚ += sum(x->-log1p(x), σcalcₚ)) 

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - ll) / simannealT(n, T0annealing, λannealing) 

                # Update jumping distribution based on size of current accepted p_move
                if dynamicjumping && r < p_move
                    if path.agepointsₚ[k] != path.agepoints[k]
                        path.σⱼtₚ[k] = ℯ * abs(path.agepointsₚ[k] - path.agepoints[k])
                    end
                    if path.Tpointsₚ[k] != path.Tpoints[k]
                        path.σⱼTₚ[k] = ℯ * abs(path.Tpointsₚ[k] - path.Tpoints[k])
                    end
                end

                # Update the currently accepted proposal
                acceptproposal!(path)
                adm = admₚ
                zdm = zdmₚ
                ll = llₚ
                npoints = npointsₚ
                copyto!(μcalc, μcalcₚ)
                copyto!(σcalc, σcalcₚ)
            end

            # Update progress meter every `progress_interval` steps
            (mod(n, progress_interval) == 0) && update!(bprogress, n)
        end
        finish!(bprogress)

        # distributions to populate
        tpointdist = fill(T(NaN), totalpoints, nsteps)
        Tpointdist = fill(T(NaN), totalpoints, nsteps)
        resultdist = fill(T(NaN), length(data), nsteps)
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
                enoughpoints = min(pointsininterval(path.agepoints, npoints, detail.agemin, detail.agemax, dt), detail.minpoints)::Int
            end
            @label crestart

            # Copy proposal from last accepted solution
            resetproposal!(path)
            admₚ = adm
            zdmₚ = zdm
            npointsₚ = npoints
            copyto!(σcalcₚ, σcalc)

            # Randomly choose an option and point (if applicable) to adjust
            r = rand()
            k = rand(Base.OneTo(npoints))

            # Adjust the proposal
            if r < p_move
                # Move one t-T point
                movepoint!(path, k)

            elseif (r < p_move+p_birth) && (npoints < maxpoints)
                # Birth: add a new model point
                k = npointsₚ = npoints + 1
                addpoint!(path, k)

            elseif (r < p_move+p_birth+p_death) && (r >= p_move+p_birth) && (npoints > max(minpoints, detail.minpoints))
                # Death: remove a model point
                npointsₚ = npoints - 1
                replacepoint!(path, k, npoints)

            elseif (r < p_move+p_birth+p_death+p_bounds)
                # Move the temperatures of the starting and ending boundaries
                # If there's an imposed unconformity or other t-T constraint, adjust within bounds
                movebounds!(path)

                # Move uncertainties
                dynamicsigma && movesigma!(σcalcₚ, data)

            elseif (r < p_move+p_birth+p_death+p_bounds+p_kinetics)
                # Adjust kinetic parameters, one at a time
                haszhe && (zdmₚ = movekinetics(zdm))
                hasahe && (admₚ = movekinetics(adm))

            end

            # Recalculate interpolated proposed t-T path
            collectproposal!(path, npointsₚ)
            
            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (pointsininterval(path.agepointsₚ, npointsₚ, detail.agemin, detail.agemax, dt) < enoughpoints) && @goto crestart
            end

            # Calculate model ages for each grain, log likelihood of proposal
            llₚ = model!(μcalcₚ, σcalcₚ, data, path.Tsteps, zdmₚ, admₚ, zftm, mftm, aftm; trackhist)
            llₚ += diff_ll(path.Tsteps, dTmax, dTmax_sigma)
            llₚ += kinetic_ll(admₚ, adm₀) + kinetic_ll(zdmₚ, zdm₀)
            simplified && (llₚ += -log(npointsₚ))
            dynamicsigma && (llₚ += sum(x->-log1p(x), σcalcₚ)) 

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - ll)

                # Update jumping distribution based on size of current accepted p_move
                if dynamicjumping && r < p_move
                    if path.agepointsₚ[k] != path.agepoints[k]
                        path.σⱼtₚ[k] = ℯ * abs(path.agepointsₚ[k] - path.agepoints[k])
                    end
                    if path.Tpointsₚ[k] != path.Tpoints[k]
                        path.σⱼTₚ[k] = ℯ * abs(path.Tpointsₚ[k] - path.Tpoints[k])
                    end
                end

                # Update the currently accepted proposal
                acceptproposal!(path)
                adm = admₚ
                zdm = zdmₚ
                ll = llₚ
                npoints = npointsₚ
                copyto!(μcalc, μcalcₚ)
                copyto!(σcalc, σcalcₚ)

                # Not critical to the function of the MCMC loop, but critical for recording stationary distribution!
                acceptancedist[n] = true
            end

            # Record results for analysis and troubleshooting
            lldist[n] = ll
            ndist[n] = npoints # distribution of # of points
            σⱼtdist[n] = path.σⱼt[k]
            σⱼTdist[n] = path.σⱼT[k]
            resultdist[:,n] .= μcalc # distribution of He ages
            admdist[n] = adm
            zdmdist[n] = zdm

            # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
            collectto!(view(tpointdist, :, n), view(path.agepoints, Base.OneTo(npoints)), path.boundary.agepoints, path.constraint.agepoints)
            collectto!(view(Tpointdist, :, n), view(path.Tpoints, Base.OneTo(npoints)), path.boundary.Tpoints, path.constraint.Tpoints)

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


    """
    ```julia
    image_from_paths(tT::TTResult; 
        xresolution::Int=1800, 
        yresolution::Int=1200, 
        xrange = nanextrema(tT.tpointdist), 
        yrange = nanextrema(tT.Tpointdist), 
    )
    ```
    Produce a 2d image (histogram) of path densities given a `TTResult`
    
    ```julia
    julia> imgcounts, xq, yq = image_from_paths(tT; xrange=boundary.agepoints, yrange=boundary.T₀)
    ```
    """
    function image_from_paths(tT::TTResult;
            xresolution::Int = 1800, 
            yresolution::Int = 1200, 
            xrange = nanextrema(tT.tpointdist), 
            yrange = nanextrema(tT.Tpointdist), 
        )
        image_from_paths(tT.tpointdist, tT.Tpointdist; xresolution, yresolution, xrange, yrange)
    end
    function image_from_paths!(tT::TTResult;
            xresolution::Int = 1800, 
            yresolution::Int = 1200, 
            xrange = nanextrema(tT.tpointdist), 
            yrange = nanextrema(tT.Tpointdist), 
        )
        image_from_paths!(tT.tpointdist, tT.Tpointdist; xresolution, yresolution, xrange, yrange)
    end
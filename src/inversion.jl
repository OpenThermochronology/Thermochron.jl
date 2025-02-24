    """
    ```julia
    MCMC(chrons::Vector{<:Chronometer}, model::NamedTuple, npoints::Int, path.agepoints::Vector, path.Tpoints::Vector, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])
    ```
    Markov chain Monte Carlo time-Temperature inversion of the thermochronometric chrons 
    specified as a vector `chrons` of `Chronometer` objects (`ZirconHe`, `ApatiteHe`, 
    `ApatiteFT`, etc.) and model parameters specified by the named tuple `model`, 
    with variable diffusion kinetics.

    Returns a `TTResult` object containing posterior time-temperature paths,

    See also `MCMC_varkinetics` for a variant with variable diffusion kinetics.

    ## Examples
    ```julia
    tT = MCMC(chrons::NamedTuple, model::NamedTuple, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])
    ```
    """
    function MCMC(dataset::NamedTuple, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where {T <: AbstractFloat}
        chrons, damodels = chronometers(T, dataset, model)
        return MCMC(chrons, damodels, model, boundary, constraint, detail)
    end
    function MCMC(chrons::Vector{<:ChronometerUnion{T}}, damodels::Vector{<:ModelUnion{T}}, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where T <: AbstractFloat
        # Process inputs
        burnin = (haskey(model, :burnin) ? model.burnin : 5*10^5)::Int
        nsteps = (haskey(model, :nsteps) ? model.nsteps : 10^6)::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        maxpoints = (haskey(model, :maxpoints) ? model.maxpoints : 50)::Int
        npoints = (haskey(model, :npoints) ? model.npoints : minpoints)::Int
        npoints = max(npoints, detail.minpoints+1)
        totalpoints = maxpoints + boundary.npoints + constraint.npoints::Int
        rescale = (haskey(model, :rescale) ? model.rescale : false)::Bool
        trackhist = (haskey(model, :trackhist) ? model.trackhist : false)::Bool
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicsigma = (haskey(model, :dynamicsigma) ? model.dynamicsigma : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : dTmax/4)::T
        agesteps = floatrange(model.agesteps)
        tsteps = floatrange(model.tsteps)
        @assert tsteps == reverse(agesteps)
        @assert issorted(tsteps)
        σmodel = T(haskey(model, :σmodel) ? model.σmodel : 1)::T
        σcalc = (haskey(model, :σcalc) ? model.σcalc : fill(σmodel, length(chrons)))::Vector{T}
        μcalc = zeros(T, length(chrons))::Vector{T}
        @assert eachindex(σcalc) == eachindex(μcalc) == eachindex(chrons)
        T0annealing = T(haskey(model, :T0annealing) ? model.T0annealing : 1)::T
        λannealing = T(haskey(model, :λannealing) ? model.λannealing : 7/burnin)::T

        # See what minerals we have
        (haszhe = any(x->isa(x, ZirconHe), chrons)) && @info "Inverting for He ages of $(count(x->isa(x, ZirconHe), chrons)) zircons"
        (hasahe = any(x->isa(x, ApatiteHe), chrons)) && @info "Inverting for He ages of $(count(x->isa(x, ApatiteHe), chrons)) apatites"
        (any(x->isa(x, SphericalHe), chrons)) && @info "Inverting for He ages of $(count(x->isa(x, SphericalHe), chrons)) generic spherical He chronometers"
        (any(x->isa(x, PlanarHe), chrons)) && @info "Inverting for He ages of $(count(x->isa(x, PlanarHe), chrons)) generic planar slab He chronometers"
        (any(x->isa(x, ZirconFT), chrons)) && @info "Inverting for fission track ages of $(count(x->isa(x, ZirconFT), chrons)) zircons"
        (any(x->isa(x, MonaziteFT), chrons)) && @info "Inverting for fission track ages of $(count(x->isa(x, MonaziteFT), chrons)) monazites"
        (any(x->isa(x, ApatiteFT), chrons)) && @info "Inverting for fission track ages of $(count(x->isa(x, ApatiteFT), chrons)) apatites"
        (any(x->isa(x, ApatiteTrackLength), chrons)) && @info "Inverting for track lengths of $(count(x->isa(x, ApatiteTrackLength), chrons)) apatite fission tracks"
        (any(x->isa(x, SphericalAr), chrons)) && @info "Inverting for Ar ages of $(count(x->isa(x, SphericalAr), chrons)) generic spherical Ar chronometers"
        (any(x->isa(x, PlanarAr), chrons)) && @info "Inverting for Ar ages of $(count(x->isa(x, PlanarAr), chrons)) generic planar slab Ar chronometers"
        (any(x->isa(x, MultipleDomain), chrons)) && @info "Inverting for ages of $(count(x->isa(x, MultipleDomain), chrons)) multiple domain diffusion chronometers"

        # Struct to hold t-T path proposals and related variables
        path = TtPath(agesteps, constraint, boundary, detail, maxpoints)

        # Variables to hold proposals
        npointsₚ = npoints
        μcalcₚ = copy(μcalc)::Vector{T}
        σcalcₚ = copy(σcalc)::Vector{T}
        
        # Initial propopsal
        initialproposal!(path, npoints) 

        # Log-likelihood for initial proposal
        ll = llₚ = model!(μcalc, σcalc, chrons, damodels, path.Tsteps; trackhist, rescale) + diff_ll(path.Tsteps, dTmax, dTmax_sigma) + 
            (simplified ? -log(npoints) : zero(T))  + (dynamicsigma ? sum(x->-log1p(x), σcalc) : zero(T)) 

        # Proposal probabilities (must sum to 1)
        p_move = 0.64
        p_birth = 0.15 # Must equal p_death
        p_death = 0.15 # Must equal p_birth
        p_bounds = 0.06

        bprogress = Progress(burnin, dt=1, desc="MCMC burn-in ($(burnin) steps)")
        progress_interval = ceil(Int,sqrt(burnin))
        for n = 1:burnin
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
                dynamicsigma && movesigma!(σcalcₚ, chrons)

            end

            # Recalculate interpolated proposed t-T path
            collectproposal!(path, npointsₚ)

            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (proposeddetail(path, npointsₚ, detail) < mindetail(path, npoints, detail)) && @goto brestart
            end

            # Calculate model ages for each grain, log likelihood of proposal
            llₚ = model!(μcalcₚ, σcalcₚ, chrons, damodels, path.Tsteps; trackhist, rescale)
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
        resultdist = fill(T(NaN), length(chrons), nsteps)
        σⱼtdist = zeros(T, nsteps)
        σⱼTdist = zeros(T, nsteps)
        lldist = zeros(T, nsteps)
        ndist = zeros(Int, nsteps)
        acceptancedist = falses(nsteps)

        progress = Progress(nsteps, dt=1, desc="MCMC collection ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps
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
                dynamicsigma && movesigma!(σcalcₚ, chrons)

            end

            # Recalculate interpolated proposed t-T path
            collectproposal!(path, npointsₚ)
            
            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (proposeddetail(path, npointsₚ, detail) < mindetail(path, npoints, detail)) && @goto crestart
            end

            # Calculate model ages for each grain, log likelihood of proposal
            llₚ = model!(μcalcₚ, σcalcₚ, chrons, damodels, path.Tsteps; trackhist, rescale)
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
    Markov chain Monte Carlo time-Temperature inversion of the thermochronometric chrons 
    specified as a vector `chrons` of `Chronometer` objects (`ZirconHe`, `ApatiteHe`, 
    `ApatiteFT`, etc.) and model parameters specified by the named tuple `model`, 
    with variable diffusion kinetics.

    Returns a `TTResult` object containing posterior time-temperature paths,
    and a `KineticResult` object containing the posterior kinetic parameters.

    See also `MCMC` for a variant with constant diffusion kinetics.

    ## Examples
    ```julia
    tT, kinetics = MCMC_varkinetics(chrons::NamedTuple, model::NamedTuple, constraint::Constraint, boundary::Boundary, [detail::DetailInterval])
    ```
    """
    function MCMC_varkinetics(dataset::NamedTuple, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where {T <: AbstractFloat}
        chrons, damodels = chronometers(T, dataset, model)
        return MCMC_varkinetics(chrons, damodels, model, boundary, constraint, detail)
    end
    function MCMC_varkinetics(chrons::Vector{<:ChronometerUnion{T}}, damodels::Vector{<:ModelUnion{T}}, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T)) where T <: AbstractFloat
        # Process inputs
        burnin = (haskey(model, :burnin) ? model.burnin : 5*10^5)::Int
        nsteps = (haskey(model, :nsteps) ? model.nsteps : 10^6)::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        maxpoints = (haskey(model, :maxpoints) ? model.maxpoints : 50)::Int
        npoints = (haskey(model, :npoints) ? model.npoints : minpoints)::Int
        npoints = max(npoints, detail.minpoints+1)
        totalpoints = maxpoints + boundary.npoints + constraint.npoints::Int
        rescale = (haskey(model, :rescale) ? model.rescale : false)::Bool
        trackhist = (haskey(model, :trackhist) ? model.trackhist : false)::Bool
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicsigma = (haskey(model, :dynamicsigma) ? model.dynamicsigma : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : dTmax/4)::T
        agesteps = floatrange(model.agesteps)
        tsteps = floatrange(model.tsteps)
        @assert tsteps == reverse(agesteps)
        @assert issorted(tsteps)
        σmodel = T(haskey(model, :σmodel) ? model.σmodel : 1)::T
        σcalc = (haskey(model, :σcalc) ? model.σcalc : fill(σmodel, length(chrons)))::Vector{T}
        μcalc = zeros(T, length(chrons))::Vector{T}
        @assert eachindex(σcalc) == eachindex(μcalc) == eachindex(chrons)
        T0annealing = T(haskey(model, :T0annealing) ? model.T0annealing : 1)::T
        λannealing = T(haskey(model, :λannealing) ? model.λannealing : 7/burnin)::T

        # See what minerals we have
        (haszhe = any(x->isa(x, ZirconHe), chrons)) && @info "Inverting for He ages of $(count(x->isa(x, ZirconHe), chrons)) zircons"
        (hasahe = any(x->isa(x, ApatiteHe), chrons)) && @info "Inverting for He ages of $(count(x->isa(x, ApatiteHe), chrons)) apatites"
        (any(x->isa(x, SphericalHe), chrons)) && @info "Inverting for He ages of $(count(x->isa(x, SphericalHe), chrons)) generic spherical He chronometers"
        (any(x->isa(x, PlanarHe), chrons)) && @info "Inverting for He ages of $(count(x->isa(x, PlanarHe), chrons)) generic planar slab He chronometers"
        (any(x->isa(x, ZirconFT), chrons)) && @info "Inverting for fission track ages of $(count(x->isa(x, ZirconFT), chrons)) zircons"
        (any(x->isa(x, MonaziteFT), chrons)) && @info "Inverting for fission track ages of $(count(x->isa(x, MonaziteFT), chrons)) monazites"
        (any(x->isa(x, ApatiteFT), chrons)) && @info "Inverting for fission track ages of $(count(x->isa(x, ApatiteFT), chrons)) apatites"
        (any(x->isa(x, ApatiteTrackLength), chrons)) && @info "Inverting for track lengths of $(count(x->isa(x, ApatiteTrackLength), chrons)) apatite fission tracks"
        (any(x->isa(x, SphericalAr), chrons)) && @info "Inverting for Ar ages of $(count(x->isa(x, SphericalAr), chrons)) generic spherical Ar chronometers"
        (any(x->isa(x, PlanarAr), chrons)) && @info "Inverting for Ar ages of $(count(x->isa(x, PlanarAr), chrons)) generic planar slab Ar chronometers"
        (any(x->isa(x, MultipleDomain), chrons)) && @info "Inverting for ages of $(count(x->isa(x, MultipleDomain), chrons)) multiple domain diffusion chronometers"

        # Struct to hold t-T path proposals and related variables
        path = TtPath(agesteps, constraint, boundary, detail, maxpoints)

        # Variables to hold proposals
        npointsₚ = npoints
        μcalcₚ = copy(μcalc)::Vector{T}
        σcalcₚ = copy(σcalc)::Vector{T}
        damodelsₚ = copy(damodels)
        damodels₀ = copy(damodels)
        updatekinetics = falses(size(damodels))

        # Initial propopsal
        initialproposal!(path, npoints)

        # Log-likelihood for initial proposal
        ll = llₚ = model!(μcalc, σcalc, chrons, damodels, path.Tsteps; trackhist, rescale) + 
            diff_ll(path.Tsteps, dTmax, dTmax_sigma) + kinetic_ll(damodelsₚ, damodels₀) + 
            (simplified ? -log(npoints) : zero(T)) + (dynamicsigma ? sum(x->-log1p(x), σcalc) : zero(T)) 
        
        # Proposal probabilities (must sum to 1)
        p_move = 0.6
        p_birth = 0.15 # Must equal p_death
        p_death = 0.15 # Must equal p_birth
        p_bounds = 0.06
        p_kinetics = 0.04

        bprogress = Progress(burnin, dt=1, desc="MCMC burn-in ($(burnin) steps)")
        progress_interval = ceil(Int,sqrt(burnin))
        for n = 1:burnin
            @label brestart

            # Copy proposal from last accepted solution
            resetproposal!(path)
            copyto!(damodelsₚ, damodels)
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
                dynamicsigma && movesigma!(σcalcₚ, chrons)

            elseif (r < p_move+p_birth+p_death+p_bounds+p_kinetics)
                # Adjust kinetic parameters, one at a time
                movekinetics!(damodelsₚ, updatekinetics)

            end

            # Recalculate interpolated proposed t-T path
            collectproposal!(path, npointsₚ)

            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (proposeddetail(path, npointsₚ, detail) < mindetail(path, npoints, detail)) && @goto brestart
            end
               
            # Calculate model ages for each grain, log likelihood of proposal
            llₚ = model!(μcalcₚ, σcalcₚ, chrons, damodelsₚ, path.Tsteps; trackhist, rescale)
            llₚ += diff_ll(path.Tsteps, dTmax, dTmax_sigma)
            llₚ += kinetic_ll(damodelsₚ, damodels₀)
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
                copyto!(damodels, damodelsₚ)
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
        resultdist = fill(T(NaN), length(chrons), nsteps)
        dmdist = similar(damodels, length(unique(damodels)), nsteps)
        σⱼtdist = zeros(T, nsteps)
        σⱼTdist = zeros(T, nsteps)
        lldist = zeros(T, nsteps)
        ndist = zeros(Int, nsteps)
        acceptancedist = falses(nsteps)

        progress = Progress(nsteps, dt=1, desc="MCMC collection ($(nsteps) steps):")
        progress_interval = ceil(Int,sqrt(nsteps))
        for n = 1:nsteps
            @label crestart

            # Copy proposal from last accepted solution
            resetproposal!(path)
            copyto!(damodelsₚ, damodels)
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
                dynamicsigma && movesigma!(σcalcₚ, chrons)

            elseif (r < p_move+p_birth+p_death+p_bounds+p_kinetics)
                # Adjust kinetic parameters, one at a time
                movekinetics!(damodelsₚ, updatekinetics)

            end

            # Recalculate interpolated proposed t-T path
            collectproposal!(path, npointsₚ)
            
            # Ensure we have enough points in the "detail" interval, if any
            if detail.minpoints > 0
                (proposeddetail(path, npointsₚ, detail) < mindetail(path, npoints, detail)) && @goto crestart
            end

            # Calculate model ages for each grain, log likelihood of proposal
            llₚ = model!(μcalcₚ, σcalcₚ, chrons, damodelsₚ, path.Tsteps; trackhist, rescale)
            llₚ += diff_ll(path.Tsteps, dTmax, dTmax_sigma)
            llₚ += kinetic_ll(damodelsₚ, damodels₀)
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
                copyto!(damodels, damodelsₚ)
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
            copyunique!(@views(dmdist[:,n]), damodels)

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
            dmdist,
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
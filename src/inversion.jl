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
    function MCMC(dataset::NamedTuple, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T); kwargs...) where {T <: AbstractFloat}
        chrons, damodels = chronometers(T, dataset, model)
        MCMC(chrons, damodels, model, boundary, constraint, detail; kwargs...)
    end
    function MCMC(chrons::Vector{<:Chronometer{T}}, damodels::Vector{<:Model{T}}, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T); liveplot::Bool=false) where T <: AbstractFloat
        # Process inputs
        burnin = (haskey(model, :burnin) ? model.burnin : 5*10^5)::Int
        nsteps = (haskey(model, :nsteps) ? model.nsteps : 10^6)::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        maxpoints = (haskey(model, :maxpoints) ? model.maxpoints : 50)::Int
        npoints = (haskey(model, :npoints) ? model.npoints : minpoints)::Int
        npoints = max(npoints, detail.minpoints+1)
        totalpoints = maxpoints + boundary.npoints + constraint.npoints::Int
        rescalemdd = (haskey(model, :rescalemdd) ? model.rescalemdd : true)::Bool
        rescale = (haskey(model, :rescale) ? model.rescale : false)::Bool
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicsigma = (haskey(model, :dynamicsigma) ? model.dynamicsigma : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : dTmax/4)::T
        agesteps = applyeltype(T, model.agesteps)
        @assert issorted(agesteps, lt=<=, rev=true) "`agesteps` must be in strictly decreasing order"
        @assert last(agesteps) >= 0 "all `agesteps` must be positive"
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
        (any(x->isa(x, ApatiteTrackLengthOriented), chrons)) && @info "Inverting for track lengths of $(count(x->isa(x, ApatiteTrackLengthOriented), chrons)) oriented apatite fission tracks"
        (any(x->isa(x, ApatiteTrackLength), chrons)) && @info "Inverting for track lengths of $(count(x->isa(x, ApatiteTrackLength), chrons)) apatite fission tracks"
        (any(x->isa(x, ZirconTrackLength), chrons)) && @info "Inverting for track lengths of $(count(x->isa(x, ZirconTrackLength), chrons)) zircon fission tracks"
        (any(x->isa(x, MonaziteTrackLength), chrons)) && @info "Inverting for track lengths of $(count(x->isa(x, MonaziteTrackLength), chrons)) monazite fission tracks"
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
        ll = llₚ = model!(μcalc, σcalc, chrons, damodels, path.Tsteps; rescalemdd, rescale) + diff_ll(path.Tsteps, dTmax, dTmax_sigma) + 
            (simplified ? -log(npoints) : zero(T))  + (dynamicsigma ? sum(x->-log1p(x), σcalc) : zero(T)) 

        # Proposal probabilities (must sum to 1)
        p_move = 0.64
        p_birth = 0.15 # Must equal p_death
        p_death = 0.15 # Must equal p_birth
        p_bounds = 0.06
        @assert p_move + p_birth + p_death + p_bounds ≈ 1
        @assert p_birth == p_death

        # Prepare and run burnin
        bprogress = Progress(burnin, desc="MCMC burn-in ($(burnin) steps)")
        progress_interval = ceil(Int,sqrt(burnin))
        if liveplot # Optionally prepare to plot t-T paths live
            imgcounts = zeros(200, 300)
            h = plot(framestyle=:box,
                xlabel="Time [Ma]", 
                ylabel="Temperature [°C]", 
                title="Paths not recorded during burnin",
                colorbar_title="Number of paths",
            )
            display(h)
            tpointbuffer = fill(T(NaN), totalpoints, progress_interval)
            Tpointbuffer = fill(T(NaN), totalpoints, progress_interval)
        end
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
            llₚ = model!(μcalcₚ, σcalcₚ, chrons, damodels, path.Tsteps; rescalemdd, rescale)
            llₚ += diff_ll(path.Tsteps, dTmax, dTmax_sigma)
            simplified && (llₚ += -log(npointsₚ))
            dynamicsigma && (llₚ += sum(x->-log1p(x), σcalcₚ)) 

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - ll) / simannealT(n, T0annealing, λannealing)
                # Update the currently accepted proposal
                dynamicjumping && r < p_move && updatejumping!(path, k)
                acceptproposal!(path)
                ll = llₚ
                npoints = npointsₚ
                copyto!(μcalc, μcalcₚ)
                copyto!(σcalc, σcalcₚ)
            end

            # Update progress meter every `progress_interval` steps
            (mod(n, progress_interval) == 0) && update!(bprogress, n)

            # Optionally update live t-T plot
            if liveplot
                tpbv = fill!(view(tpointbuffer, :, 1+mod(n, progress_interval)), NaN)
                Tpbv = fill!(view(Tpointbuffer, :, 1+mod(n, progress_interval)), NaN)
                collectto!(tpbv, view(path.agepoints, Base.OneTo(npoints)), path.boundary.agepoints, path.constraint.agepoints)
                collectto!(Tpbv, view(path.Tpoints, Base.OneTo(npoints)), path.boundary.Tpoints, path.constraint.Tpoints)
                if mod(n, progress_interval) == 0
                    (A, xc, yc) = image_from_paths!(imgcounts, tpointbuffer, Tpointbuffer; xrange=boundary.agepoints, yrange=boundary.T₀)
                    heatmap!(h, xc, yc, A, colormap=:viridis, zlims=(0, nanpctile(A,95)), title="Burn-in: $n of $nsteps steps")
                    display(h)
                end
            end
        end
        finish!(bprogress)
        liveplot && fill!(imgcounts, 0)

        # distributions to populate
        tpointdist = fill(T(NaN), totalpoints, nsteps)
        Tpointdist = fill(T(NaN), totalpoints, nsteps)
        resultdist = fill(T(NaN), length(chrons), nsteps)
        σⱼtdist = zeros(T, nsteps)
        σⱼTdist = zeros(T, nsteps)
        lldist = zeros(T, nsteps)
        ndist = zeros(Int, nsteps)
        acceptancedist = falses(nsteps)

        progress = Progress(nsteps, desc="MCMC collection ($(nsteps) steps):")
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
            llₚ = model!(μcalcₚ, σcalcₚ, chrons, damodels, path.Tsteps; rescalemdd, rescale)
            llₚ += diff_ll(path.Tsteps, dTmax, dTmax_sigma)
            simplified && (llₚ += -log(npointsₚ))
            dynamicsigma && (llₚ += sum(x->-log1p(x), σcalcₚ)) 

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - ll)
                # Update the currently accepted proposal
                dynamicjumping && r < p_move && updatejumping!(path, k)
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

            # Optionally update t-T plot
            if liveplot && mod(n, progress_interval) == 0
                new = (n-progress_interval+1):n
                (A, xc, yc) = image_from_paths!(imgcounts, tpointdist[:,new], Tpointdist[:,new]; xrange=boundary.agepoints, yrange=boundary.T₀)
                heatmap!(h, xc, yc, A, colormap=:viridis, zlims=(0, nanpctile(A,95)), title="Collection: $n of $nsteps steps")
                display(h)
            end
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
    function MCMC_varkinetics(dataset::NamedTuple, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T); kwargs...) where {T <: AbstractFloat}
        chrons, damodels = chronometers(T, dataset, model)
        MCMC_varkinetics(chrons, damodels, model, boundary, constraint, detail; kwargs...)
    end
    function MCMC_varkinetics(chrons::Vector{<:Chronometer{T}}, damodels::Vector{<:Model{T}}, model::NamedTuple, boundary::Boundary{T}, constraint::Constraint{T}=Constraint(T), detail::DetailInterval{T}=DetailInterval(T); liveplot::Bool=false) where T <: AbstractFloat
        # Process inputs
        burnin = (haskey(model, :burnin) ? model.burnin : 5*10^5)::Int
        nsteps = (haskey(model, :nsteps) ? model.nsteps : 10^6)::Int
        minpoints = (haskey(model, :minpoints) ? model.minpoints : 1)::Int
        maxpoints = (haskey(model, :maxpoints) ? model.maxpoints : 50)::Int
        npoints = (haskey(model, :npoints) ? model.npoints : minpoints)::Int
        npoints = max(npoints, detail.minpoints+1)
        totalpoints = maxpoints + boundary.npoints + constraint.npoints::Int
        rescalemdd = (haskey(model, :rescalemdd) ? model.rescalemdd : true)::Bool
        rescale = (haskey(model, :rescale) ? model.rescale : false)::Bool
        simplified = (haskey(model, :simplified) ? model.simplified : false)::Bool
        dynamicsigma = (haskey(model, :dynamicsigma) ? model.dynamicsigma : false)::Bool
        dynamicjumping = (haskey(model, :dynamicjumping) ? model.dynamicjumping : false)::Bool
        dTmax = T(haskey(model, :dTmax) ? model.dTmax : 10)::T
        dTmax_sigma = T(haskey(model, :dTmax_sigma) ? model.dTmax_sigma : dTmax/4)::T
        agesteps = applyeltype(T, model.agesteps)
        @assert issorted(agesteps, lt=<=, rev=true) "`agesteps` must be in strictly decreasing order"
        @assert last(agesteps) >= 0 "all `agesteps` must be positive"
        σmodel = T(haskey(model, :σmodel) ? model.σmodel : 1)::T
        σcalc = (haskey(model, :σcalc) ? model.σcalc : fill(σmodel, length(chrons)))::Vector{T}
        μcalc = zeros(T, length(chrons))::Vector{T}
        @assert eachindex(σcalc) == eachindex(μcalc) == eachindex(chrons)
        T0annealing = T(haskey(model, :T0annealing) ? model.T0annealing : 1)::T
        λannealing = T(haskey(model, :λannealing) ? model.λannealing : 7/burnin)::T
        redegasparent = true

        # See what minerals we have
        (haszhe = any(x->isa(x, ZirconHe), chrons)) && @info "Inverting for He ages of $(count(x->isa(x, ZirconHe), chrons)) zircons"
        (hasahe = any(x->isa(x, ApatiteHe), chrons)) && @info "Inverting for He ages of $(count(x->isa(x, ApatiteHe), chrons)) apatites"
        (any(x->isa(x, SphericalHe), chrons)) && @info "Inverting for He ages of $(count(x->isa(x, SphericalHe), chrons)) generic spherical He chronometers"
        (any(x->isa(x, PlanarHe), chrons)) && @info "Inverting for He ages of $(count(x->isa(x, PlanarHe), chrons)) generic planar slab He chronometers"
        (any(x->isa(x, ZirconFT), chrons)) && @info "Inverting for fission track ages of $(count(x->isa(x, ZirconFT), chrons)) zircons"
        (any(x->isa(x, MonaziteFT), chrons)) && @info "Inverting for fission track ages of $(count(x->isa(x, MonaziteFT), chrons)) monazites"
        (any(x->isa(x, ApatiteFT), chrons)) && @info "Inverting for fission track ages of $(count(x->isa(x, ApatiteFT), chrons)) apatites"
        (any(x->isa(x, ApatiteTrackLengthOriented), chrons)) && @info "Inverting for track lengths of $(count(x->isa(x, ApatiteTrackLengthOriented), chrons)) oriented apatite fission tracks"
        (any(x->isa(x, ApatiteTrackLength), chrons)) && @info "Inverting for track lengths of $(count(x->isa(x, ApatiteTrackLength), chrons)) apatite fission tracks"
        (any(x->isa(x, ZirconTrackLength), chrons)) && @info "Inverting for track lengths of $(count(x->isa(x, ZirconTrackLength), chrons)) zircon fission tracks"
        (any(x->isa(x, MonaziteTrackLength), chrons)) && @info "Inverting for track lengths of $(count(x->isa(x, MonaziteTrackLength), chrons)) monazite fission tracks"
        (any(x->isa(x, SphericalAr), chrons)) && @info "Inverting for Ar ages of $(count(x->isa(x, SphericalAr), chrons)) generic spherical Ar chronometers"
        (any(x->isa(x, PlanarAr), chrons)) && @info "Inverting for Ar ages of $(count(x->isa(x, PlanarAr), chrons)) generic planar slab Ar chronometers"
        (any(x->isa(x, MultipleDomain), chrons)) && @info "Inverting for ages of $(count(x->isa(x, MultipleDomain), chrons)) multiple domain diffusion chronometers"

        # Struct to hold t-T path proposals and related variables
        path = TtPath(agesteps, constraint, boundary, detail, maxpoints)

        # Variables to hold proposals
        npointsₚ = npoints
        μcalcₚ = copy(μcalc)::Vector{T}
        σcalcₚ = copy(σcalc)::Vector{T}
        damodels₀ = damodels
        damodels = copy(damodels)
        damodelsₚ = copy(damodels)
        updatekinetics = falses(size(damodels))

        # Initial propopsal
        initialproposal!(path, npoints)

        # Log-likelihood for initial proposal
        ll = llₚ = model!(μcalc, σcalc, chrons, damodels, path.Tsteps; rescalemdd, rescale, redegasparent) + 
            diff_ll(path.Tsteps, dTmax, dTmax_sigma) + kinetic_ll(damodelsₚ, damodels₀) + 
            (simplified ? -log(npoints) : zero(T)) + (dynamicsigma ? sum(x->-log1p(x), σcalc) : zero(T)) 
        
        # Proposal probabilities (must sum to 1)
        p_move = 0.6
        p_birth = 0.14 # Must equal p_death
        p_death = 0.14 # Must equal p_birth
        p_bounds = 0.06
        p_kinetics = 0.06
        @assert p_move + p_birth + p_death + p_bounds + p_kinetics ≈ 1
        @assert p_birth == p_death

        # Prepare and run burnin
        bprogress = Progress(burnin, desc="MCMC burn-in ($(burnin) steps)")
        progress_interval = ceil(Int,sqrt(burnin))
        if liveplot # Optionally prepare to plot t-T paths
            imgcounts = zeros(200, 300)
            h = plot(framestyle=:box,
                xlabel="Time [Ma]", 
                ylabel="Temperature [°C]", 
                colorbar_title="Number of paths",
            )
            display(h)
            tpointbuffer = fill(T(NaN), totalpoints, progress_interval)
            Tpointbuffer = fill(T(NaN), totalpoints, progress_interval)
        end
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
            llₚ = model!(μcalcₚ, σcalcₚ, chrons, damodelsₚ, path.Tsteps; rescalemdd, rescale, redegasparent)
            llₚ += diff_ll(path.Tsteps, dTmax, dTmax_sigma)
            llₚ += kinetic_ll(damodelsₚ, damodels₀)
            simplified && (llₚ += -log(npointsₚ))
            dynamicsigma && (llₚ += sum(x->-log1p(x), σcalcₚ)) 

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - ll) / simannealT(n, T0annealing, λannealing) 
                # Update the currently accepted proposal
                dynamicjumping && r < p_move && updatejumping!(path, k)
                acceptproposal!(path)
                copyto!(damodels, damodelsₚ)
                ll = llₚ
                npoints = npointsₚ
                copyto!(μcalc, μcalcₚ)
                copyto!(σcalc, σcalcₚ)
            end

            # Update progress meter every `progress_interval` steps
            (mod(n, progress_interval) == 0) && update!(bprogress, n)

            # Optionally update live t-T plot
            if liveplot
                tpbv = fill!(view(tpointbuffer, :, 1+mod(n, progress_interval)), NaN)
                Tpbv = fill!(view(Tpointbuffer, :, 1+mod(n, progress_interval)), NaN)
                collectto!(tpbv, view(path.agepoints, Base.OneTo(npoints)), path.boundary.agepoints, path.constraint.agepoints)
                collectto!(Tpbv, view(path.Tpoints, Base.OneTo(npoints)), path.boundary.Tpoints, path.constraint.Tpoints)
                if mod(n, progress_interval) == 0
                    (A, xc, yc) = image_from_paths!(imgcounts, tpointbuffer, Tpointbuffer; xrange=boundary.agepoints, yrange=boundary.T₀)
                    heatmap!(h, xc, yc, A, colormap=:viridis, zlims=(0, nanpctile(A,95)), title="Burn-in: $n of $nsteps steps")
                    display(h)
                end
            end
        end
        finish!(bprogress)
        liveplot && fill!(imgcounts, 0)

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

        progress = Progress(nsteps, desc="MCMC collection ($(nsteps) steps):")
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
            llₚ = model!(μcalcₚ, σcalcₚ, chrons, damodelsₚ, path.Tsteps; rescalemdd, rescale, redegasparent)
            llₚ += diff_ll(path.Tsteps, dTmax, dTmax_sigma)
            llₚ += kinetic_ll(damodelsₚ, damodels₀)
            simplified && (llₚ += -log(npointsₚ))
            dynamicsigma && (llₚ += sum(x->-log1p(x), σcalcₚ)) 

            # Accept or reject proposal based on likelihood
            if log(rand()) < (llₚ - ll)
                # Update the currently accepted proposal
                dynamicjumping && r < p_move && updatejumping!(path, k)
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

            # Optionally update live t-T plot
            if liveplot && mod(n, progress_interval) == 0
                new = (n-progress_interval+1):n
                (A, xc, yc) = image_from_paths!(imgcounts, tpointdist[:,new], Tpointdist[:,new]; xrange=boundary.agepoints, yrange=boundary.T₀)
                heatmap!(h, xc, yc, A, colormap=:viridis, zlims=(0, nanpctile(A,95)), title="Collection: $n of $nsteps steps")
                display(h)
            end
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
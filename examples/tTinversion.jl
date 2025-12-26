# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                              tTinversion.jl                                   #
#                                                                               #
#       A transdimensional Bayesian Markov-chain Monte Carlo code for           #
#   inverting time-temperature paths from zircon and apatite U-Th/He ages       #
#   apatite fission track ages, and apatite fission track lengths.              #
#                                                                               #
#       This file uses the MCMC, damage annealing, and Crank-Nicholson          #
#   diffusion codes provided by the Thermochron.jl package.                     #
#                                                                               #
#   © 2025 C. Brenhin Keller and Kalin McDannell                                #
#                                                                               #
#       If running for the first time, you might consider instantiating the     #
#   manifest that came with this file                                           #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## ---  Load required packages

    using Thermochron, Plots
    using LinearAlgebra
    # Diminishing returns with more than ~2 threads
    BLAS.get_num_threads() > 2 && BLAS.set_num_threads(2)

    # Make sure we're running in the directory where the script is located
    cd(@__DIR__)

    # # # # # # # # # # Choice of regional thermochron data # # # # # # # # # #

    # Literature samples from McDannell et al. 2022 (doi: 10.1130/G50315.1), Manitoba
    # (12 ZirconHe, 5 ApatiteHe, 47 ApatiteFT, 269 ApatiteTrackLengthOriented)
    name = "Manitoba"
    datafile = joinpath("exampledata", "manitoba.csv")

    # # Literature samples from Guenthner 2013 (doi: 10.2475/03.2013.01), Minnesota
    # # (23 ZirconHe, 11 ApatiteHe)
    # name = "Minnesota"
    # datafile = joinpath("exampledata", "minnesota.csv")

    # # Literature samples from Valla et al. 2011 (doi: 10.1038/NGEO1242), Visp, Switzerland
    # # (23 ApatiteHe, 4 He-4/He-3)
    # name = "VIS"
    # datafile = joinpath("exampledata", "vis.csv")

    # # Literature sample from McDannell et al. 2018 (doi: 10.1016/j.epsl.2018.03.012), Quebec
    # # (1 MDD)
    # name = "OL13"
    # datafile = joinpath("exampledata", "ol13.csv")

    ds = importdataset(datafile, ',', importas=:Tuple)

## --- Prepare problem

    # Option A: Linear timesteps
    dt = 8 # [Ma] Average model timestep (smaller is slower but more accurate)
    tinit = ceil(maximum(ds.crystallization_age_Ma)/dt) * dt # [Ma] Model start time
    agesteps = cntr(tinit:-dt:0)

    # # Option B: Logarithmic timesteps
    # resolution = 300 # Approximate number of model timesteps (more is slower but more accurate)
    # dt = maximum(ds.crystallization_age_Ma)/resolution # [Ma] Average model timestep
    # tinit = maximum(ds.crystallization_age_Ma) # [Ma] Model start time
    # agesteps = cntr(logrange(tinit+2dt, 2dt, length=Int(tinit÷dt+1)) .- 2dt)

    # Check the validity of the requested time discretization (correct order, etc.)
    agesteps, tsteps = checktimediscretization(Float64, agesteps)
    @info """Preparing model with $(length(agesteps)) $(allequal(diff(agesteps)) ? "linear" : "nonlinear") timesteps with bin centers between $(round(first(agesteps),digits=3)) and $(round(last(agesteps), digits=3)) Ma"""

    params = (
        nsteps = 400000,                # [n] How many steps of the Markov chain should we run?
        burnin = 100000,                # [n] How long should we wait for MC to converge (become stationary)
        dr = 1.0,                       # [μm] Radius step size
        dTmax = 25.0,                   # [Ma/dt] Maximum reheating/burial per model timestep. If too high, may cause numerical problems in Crank-Nicholson solve
        Tinit = 400.0,                  # [C] Initial model temperature (i.e., crystallization temperature)
        ΔTinit = -100.0,                # [C] Tinit can vary from Tinit to Tinit+ΔTinit
        Tnow = 0.0,                     # [C] Current surface temperature
        ΔTnow = 20.0,                   # [C] Tnow may vary from Tnow to Tnow+ΔTnow
        dt = dt,                        # [Ma] Model timestep
        tinit = tinit,                  # [Ma] Model start time
        tnow = 0.0,                     # [Ma] Model end time (today)
        tsteps = tsteps,                # [Ma] Forward time discretiziation (from start of model)
        agesteps = agesteps,            # [Ma] Age discretization (relative to the present)
        minpoints = 15,                 # [n] Minimum allowed number of model t-T points (nodes)
        maxpoints = 50,                 # [n] Maximum allowed number of model t-T points (nodes)
        dynamicjumping = true,          # Update the t and T jumping (proposal) distributions based on previously accepted jumps
        stepwisetracerfraction = false, # Calculate degassing ll based on stepwise fraction degassed (rather than cumulative)
        # Damage and annealing models for diffusivity (specify custom kinetics if desired)
        adm = RDAAM(),                  # Flowers et al. 2009 (doi: 10.1016/j.gca.2009.01.015) apatite diffusivity model
        zdm = ZRDAAM(),                 # Guenthner et al. 2013 (doi: 10.2475/03.2013.01) zircon diffusivity model
        uaftm = Ketcham1999FC(:unoriented), # Ketcham et al. 1999 (doi: 10.2138/am-1999-0903) unoriented apatite fission track model  
        aftm = Ketcham1999FC(),         # Ketcham et al. 1999 (doi: 10.2138/am-1999-0903) oriented apatite fission track model
        zftm = Yamada2007PC(),          # Yamada et al. 2007 (doi: 10.1016/j.chemgeo.2006.09.002) zircon fission track model
        mftm = Jones2021FA(),           # Re-fit from Jones et al. 2021 (doi: 10.5194/gchron-3-89-2021) monazite fission track 
        # Optional simulated annealing during burnin, wherein p_accept = max(exp((llₚ-ll)/T), 1)
        # T = T0annealing * exp(-λannealing * n) + 1 at step number n of burnin
        T0annealing = 5,                # [unitless] initial annealing "temperature" (set to 0 for no simulated annealing).
    )

    # Default: no detail interval
    detail = DetailInterval()

    # # (Optional) Uncomment this section to require greater t-T node density 
    # # in some time interval (typically the youngest end of the total time 
    # # interval, where you may expect the data to have more resolving power)
    # detail = DetailInterval(
    #     agemin = 0.0, # Youngest end of detail interval
    #     agemax = tinit/3, # Oldest end of detail interval
    #     minpoints = 8, # Minimum number of points in detail interval
    # )

    # Boundary conditions (e.g. 10C at present and 650 C at the time of zircon formation).
    boundary = Boundary(
        agepoints = [0.0, params.tinit],          # [Ma] Final and initial time
        T₀ = [params.Tnow, params.Tinit],          # [C] Final and initial temperature
        ΔT = [params.ΔTnow, params.ΔTinit],        # [C] Final and initial temperature range (positive or negative)
        tboundary = :reflecting, # Reflecting time boundary conditions
        Tboundary = :reflecting, # Reflecting temperature boundary conditions
    )

    # Default: No constraints are imposed
    constraint = Constraint()

    # # (Optional) Uncomment this section if you wish to impose an unconformity 
    # # or other constraint at any point in the record.
    # constraint = Constraint(
    #     agedist = [ Normal(510,20),],  # [Ma] Age distribution
    #     Tdist =   [  Uniform(0,50),],  # [C] Temperature distribution
    # )
    # name *= "_constrained"

## --- Process data into Chronometer objects

    chrons, damodels = chronometers(ds, params)

    # Model uncertainty is not generally well known (depending on annealing parameters,
    # decay constants, diffusion parameters, etc.), but is certainly non-zero. In addition, 
    # observed thermochronometric ages often display excess dispersion beyond analytical
    # uncertainty. Here we specify a default minimum uncertainty representing the average expected
    # misfit between model and data in excess of that expected from analytical uncertainty.
    params = (;params...,
        σcalc = 0.025*Thermochron.value.(chrons),     # [Ma] model uncertainty (2.5%)
    )

## --- Age uncertainty resampling: estimate expected misfit (σcalc) from excess dispersion of the data itself

    # Empirical age uncertainty for apatite
    tap = isa.(chrons, ApatiteHe)
    if count(tap) > 1
        h = ageeuplot(chrons[tap], label="Internal uncertainty", title="apatite")
        empiricaluncertainty!(params.σcalc, chrons, ApatiteHe, sigma_eU=10, sigma_offset=10)
        σtotal = sqrt.(get_age_sigma(chrons[tap]).^2 + params.σcalc[tap].^2)
        ageeuplot!(h, chrons[tap], yerror=2*σtotal, label="Empirical uncertainty")
        savefig(h, name*"_ap-empirical.pdf")
        display(h)
    end

    # Empirical age uncertainty for zircon
    tzr = isa.(chrons, ZirconHe)
    if count(tzr) > 1
        h = ageeuplot(chrons[tzr], label="Internal uncertainty", title="zircon")
        empiricaluncertainty!(params.σcalc, chrons, ZirconHe, sigma_eU=100, sigma_offset=10)
        σtotal = sqrt.(get_age_sigma(chrons[tzr]).^2 + params.σcalc[tzr].^2)
        ageeuplot!(h, chrons[tzr], yerror=2*σtotal, label="Empirical uncertainty")
        savefig(h, name*"_zir-empirical.pdf")
        display(h)
    end

## --- Invert for maximum likelihood t-T path

    # Run Markov Chain
    # @time tT = MCMC(chrons, damodels, params, boundary, constraint, detail; liveplot=true)
    @time tT, kinetics = MCMC_varkinetics(chrons, damodels, params, boundary, constraint, detail; liveplot=true)
    @info """tT.tpointdist & tT.Tpointdist collected, size: $(size(tT.Tpointdist))
    Mean log-likelihood: $(nanmean(view(tT.lldist, params.burnin:params.nsteps)))
    Mean acceptance rate: $(nanmean(view(tT.acceptancedist, params.burnin:params.nsteps)))
    Mean npoints: $(nanmean(view(tT.ndist, params.burnin:params.nsteps)))
    Mean jₜ: $(nanmean(view(tT.jtdist,params.burnin:params.nsteps)))
    Mean jT: $(nanmean(view(tT.jTdist, params.burnin:params.nsteps)))
    """

## --- (optional) Save or load full tT results to/from file
    using JLD2, CodecZlib

    # Save tTs using JLD2 (compressed)
    @save "$name.jld2" {compress=true} tT kinetics params

    # # Load saved tTs
    # @load "$name.jld2"

## ---  Plot log likelihood distribution

    h = plot(tT.lldist, xlabel="Step number", ylabel="Log likelihood", label="", framestyle=:box)
    savefig(h, name*"_lldist.pdf")
    display(h)

## --- Plot distribution of number of model t-T points (nodes)

    h = plot(tT.ndist, label="", framestyle=:box, lc=:purple)
    plot!(xlabel="Step number", ylabel="Number of model t-T nodes", ylims=(0,last(ylims())))
    savefig(h, name*"_ndist.pdf")
    display(h)

## --- Plot moving average of acceptance distribution

    h = plot(movmean(tT.acceptancedist,100), label="", framestyle=:box, lc=:limegreen)
    plot!(xlabel="Step number", ylabel="Acceptance probability (mean of 100)", ylims=(0,1))
    savefig(h, name*"_acceptance.pdf")
    display(h)

## --- Plot calculated/observed ages as a function of grain size, colored by eU (ZirconHe, ApatiteHe)

    C = (ApatiteHe, ZirconHe,)
    mineral = ("apatite", "zircon",)
    for i in eachindex(C, mineral)
        t = isa.(chrons, C[i])
        if any(t)
            σtotal = sqrt.(get_age_sigma(chrons[t]).^2 + params.σcalc[t].^2)
            h = agesizeplot(chrons[t], yerror=2σtotal, zcolor=eU.(chrons[t]),
                colorbar_title = "eU [ppm]",
                title = "$(C[i]): Age-Size",
                label = "Data (2σ total)", 
                shape = :hexagon,
                color = :viridis,
                msc = :black,
                lc = :black,
            )
            agedist = tT.resultdist[t,:]
            m = nanmean(agedist, dim=2)
            l = nanpctile(agedist, 2.5, dim=2)
            u = nanpctile(agedist, 97.5, dim=2)
            mcolor = mineralcolors[mineral[i]]
            scatter!(h, radius.(chrons[t]), m, yerror=(m-l, u-m), zcolor=eU.(chrons[t]),
                label = "Model (95%CI)", 
                color = :viridis,
                msc = mcolor,
                lc = mcolor,
            )
            savefig(h, "$(name)_$(C[i])_Age-Size-eU.pdf")
            display(h)
        end
    end

## --- Plot calculated/observed ages as a function of eU (ZirconHe, ApatiteHe)

    C = (ApatiteHe, ZirconHe,)
    mineral = ("apatite", "zircon",)
    for i in eachindex(C, mineral)
        t = isa.(chrons, C[i])
        if any(t)
            σtotal = sqrt.(get_age_sigma(chrons[t]).^2 + params.σcalc[t].^2)
            h = ageeuplot(chrons[t], yerror=2σtotal,
                title = "$(C[i]): Age-eU",
                label = "Data (2σ total)", 
                shape = :hexagon,
                color = :black,
                lc = :black,
            )
            agedist = tT.resultdist[t,:]
            m = nanmean(agedist, dim=2)
            l = nanpctile(agedist, 2.5, dim=2)
            u = nanpctile(agedist, 97.5, dim=2)
            mcolor = mineralcolors[mineral[i]]
            scatter!(h, eU.(chrons[t]), m, yerror=(m-l, u-m), 
                label = "Model (95%CI)", 
                color = mcolor, 
                msc = mcolor, 
            )
            savefig(h, "$(name)_$(C[i])_Age-eU.pdf")
            display(h)
        end
    end

    if @isdefined kinetics
        for D in (RDAAM, ZRDAAM)
            id = findfirst(x->isa(x, D), kinetics.dmdist[:,1])
            im = findfirst(x->isa(x, D), damodels)
            if !isnothing(id) && !isnothing(im)
                hdm = plot(damodels[im], kinetics.dmdist[id,:])
                plot!(hdm[1], title="$D kinetics")
                savefig(hdm, "$(name)_$(D)_kinetics.pdf")
                display(hdm)
            end
        end
    end

## -- Plot calculated/observed ages as a function of rₘᵣ₀ (ApatiteFT)

t = isa.(chrons, ApatiteFT)
if any(t)
    rmr0 = chrons[t] .|> x->x.rmr0
    μobs = get_age(chrons[t])
    σtotal = sqrt.(get_age_sigma(chrons[t]).^2 + params.σcalc[t].^2)
    h = scatter(rmr0, μobs, yerror=2σtotal,
        title = "Apatite FT",
        xlabel = "rₘᵣ₀ [unitless]", ylabel = "Age [Ma]",
        label = "Data (2σ total)",
        legend = :topright,
        framestyle = :box,
        shape = :hexagon,
        color = :black,
    )
    agedist = tT.resultdist[t,:]
    m = nanmean(agedist, dim=2)
    l = nanpctile(agedist, 2.5, dim=2)
    u = nanpctile(agedist, 97.5, dim=2)
    mcolor = mineralcolors["apatite"]
    scatter!(h, rmr0, m, yerror=(m-l, u-m), 
        label = "Model (95%CI)", 
        color = mcolor, 
        msc = mcolor,
    )
    # Weighted Mean of observed ages
    weights = 1.0 ./ (get_age_sigma(chrons[t]).^2)
    wmean_obs = sum(μobs .* weights) / sum(weights)
    σ_obs = sqrt(1.0 / sum(weights))  # uncertainty of weighted mean
    
    # Predicted: 95% CI
    mean_pred = nanmean(agedist)
    ci_low = nanpctile(agedist, 2.5)
    ci_high = nanpctile(agedist, 97.5)  
    
    # Annotate weighted mean and predicted ages
    xl, yl = xlims(h), ylims(h)
    ann = [text("Observed: $(round(wmean_obs, digits=1)) ± $(round(2*σ_obs, digits=1)) Ma (2σ weighted mean)", 8, :left, :black),
           text("Predicted: $(round(mean_pred, digits=1)) +$(round(ci_high-mean_pred, digits=1))/-$(round(ci_high-mean_pred, digits=1)) Ma (95% CI)", 8, :left, mineralcolors["apatite"]),]
    for i in eachindex(ann)
        annotate!(h, (minimum(xl)+0.03*nanrange(xl), minimum(yl) + nanrange(yl) * (1 - 0.06 * i), ann[i]))
    end

    savefig(h, "$(name)_ApatiteFT_predicted.pdf")
    display(h)
end

## -- Plot calculated and observed ages (most other chronometers)

    C = (PlanarAr, PlanarHe, SphericalAr, SphericalHe, ZirconFT, MonaziteFT)
    mineral = ("feldspar", "hematite", "feldspar", "hematite", "zircon", "monazite")
    for i in eachindex(C, mineral)
        t = isa.(chrons, C[i])
        if any(t)
            μobs = get_age(chrons[t])
            σtotal = sqrt.(get_age_sigma(chrons[t]).^2 + params.σcalc[t].^2)
            h = plot(μobs, yerror=2σtotal,
                title = "$(C[i])",
                xlabel = "Sample number", ylabel = "Age [Ma]",
                label = "Data (2σ total)", 
                framestyle = :box,
                shape = :hexagon,
                color = :black,
            )
            agedist = tT.resultdist[t,:]
            m = nanmean(agedist, dim=2)
            l = nanpctile(agedist, 2.5, dim=2)
            u = nanpctile(agedist, 97.5, dim=2)
            scatter!(h, m, yerror=(m-l, u-m), 
                label = "Model (95%CI)", 
                color = mineralcolors[mineral[i]], 
                msc = mineralcolors[mineral[i]], 
            )
            savefig(h, "$(name)_$(C[i])_predicted.pdf")
            display(h)
        end
    end

## -- Fission track length histograms (apatite, zircon, monazite)
    # Use HypothesisTests package for testing equivalence of distributions a posteriori (K-S test)
    using HypothesisTests

    C = (ApatiteTrackLengthOriented, ApatiteTrackLength, ZirconTrackLength, MonaziteTrackLength)
    mineral = ("apatite", "apatite", "zircon", "monazite")
    for i in eachindex(C, mineral)
        t = isa.(chrons, C[i])
        if any(t)
            # Extract observed and modeled lengths
            obs_lengths = Thermochron.value.(chrons[t])
            pred_lengths = vec(tT.resultdist[t, :])
    
            mcolor = mineralcolors[mineral[i]]
            h = histogram(obs_lengths, bins=0:0.25:20, normalized=true,
                xlims = (0,20),
                title = "$(C[i])",
                xlabel = "Track length [μm]",
                ylabel = "Probability density",
                label = "Data (N=$(count(t)))", 
                legend = :topleft,
                framestyle = :box,
                color = :black,
                alpha = 0.75,
            )
            histogram!(h, pred_lengths, bins=0:0.25:20, normalized=true, 
                label = "Model",
                color = mcolor,
                fill = true,
                alpha = 0.75,
                lw = 0,
            )
            yl = ylims(h)
    
            # Test for homogeneity
            ks = ApproximateTwoSampleKSTest(obs_lengths, pred_lengths)
            GoF = 1.0 - ks.δ

            # Annotate statistics and GOF
            ann = [text("Observed:   $(round(nanmean(obs_lengths), digits=2)) ± $(round(2*nanstd(obs_lengths), digits=2)) µm (2σ)", 8, :left, :black),
                   text("Predicted:   $(round(nanmean(pred_lengths), digits=2)) ± $(round(2*nanstd(pred_lengths), digits=2)) µm (2σ)", 8, :left, mcolor),
                   text("(K-S) GoF = $(round(GoF, digits=2)), p = $(round(pvalue(ks), digits=2))", 8, :left, :black),]
            for i in eachindex(ann)
                annotate!(h, (0.7, maximum(yl) * (0.86 - 0.06 * i), ann[i]))
            end
            
            savefig(h, "$(name)_$(C[i])_predicted.pdf")
            display(h)
        end
    end

## -- Fission track length histograms on a "sample-by-sample" basis (apatite, zircon, monazite)
    # Samples defined on the basis of sample IDs specified in the "grain name" column in the input file
    # Uses HypothesisTests package for testing equivalence of distributions.

    sample_ids = chrons .|> x->x.name # Sample IDs as defined by "grain name" column of input file
    # sample_ids = chrons .|> x->x.notes # [alternatively] Sample IDs as defined by "notes" column of input file
    C = (ApatiteTrackLengthOriented, ApatiteTrackLength, ZirconTrackLength, MonaziteTrackLength)
    mineral = ("apatite", "apatite", "zircon", "monazite")
    for i in eachindex(C, mineral)
        # Filter chronometers
        t = isa.(chrons, C[i])
        for sid in unique(sample_ids[t])
            ts = t .& (sample_ids .== sid)
            any(ts) || continue

            # Extract observed and modeled lengths
            obs_lengths = Thermochron.value.(chrons[ts])
            pred_lengths = vec(tT.resultdist[ts, :])
    
            mcolor = mineralcolors[mineral[i]]
            h = histogram(obs_lengths, bins=0:0.25:20, normalized=true,
                xlims = (0,20),
                title = "$sid",
                xlabel = "Track length [μm]",
                ylabel = "Probability density",
                label = "Data (N=$(count(ts)))", 
                legend = :topleft,
                framestyle = :box,
                color = :black,
                alpha = 0.75,
            )
            histogram!(h, pred_lengths, bins=0:0.25:20, normalized=true, 
                label = "Model",
                color = mcolor,
                fill = true,
                alpha = 0.75,
                lw = 0,
            )
            yl = ylims(h)

            # Test for homogeneity
            ks = ApproximateTwoSampleKSTest(obs_lengths, pred_lengths)
            GoF = 1.0 - ks.δ

            # Annotate statistics and GOF
            ann = [text("Observed:   $(round(nanmean(obs_lengths), digits=2)) ± $(round(2*nanstd(obs_lengths), digits=2)) µm (2σ)", 8, :left, :black),
                   text("Predicted:   $(round(nanmean(pred_lengths), digits=2)) ± $(round(2*nanstd(pred_lengths), digits=2)) µm (2σ)", 8, :left, mcolor),
                   text("(K-S) GoF = $(round(GoF, digits=2)), p = $(round(pvalue(ks), digits=2))", 8, :left, :black),]
            for i in eachindex(ann)
                annotate!(h, (0.7, maximum(yl) * (0.86 - 0.06 * i), ann[i]))
            end

            savefig(h, "$(name)_$(C[i])_$(sid)_predicted.pdf")
            display(h)
        end
    end

## --- Plot calculated and observed step heating degassing curves

    nplot = 33
    plotindices = rand(eachindex(tT), nplot)
    C = (SingleDomain, MultipleDomain)
    D = (SDDiffusivity, MDDiffusivity)
    for (C,D) in zip((SingleDomain, MultipleDomain), (SDDiffusivity, MDDiffusivity))
        dis = if @isdefined kinetics
            findall(x->x isa D, first(kinetics))
        else
            findall(x->x isa D, damodels)
        end
        cis = findall(x->x isa C, chrons)
        @assert eachindex(cis) == eachindex(dis)
        for (ci, di) in zip(cis, dis)
            c = chrons[ci]
            h = plot(title = "$(c.name) $C degassing",
                xlabel = "Time [s]",
                ylabel = "Cumulative fraction released",
                framestyle = :box,
            )
            tmin = (findfirst(c.fit) > 1) ? c.tsteps_experimental[findfirst(c.fit)-1] : 0
            tmax = c.tsteps_experimental[findlast(c.fit)]
            t = tmin .<= c.tsteps_degassing .<= tmax
            for pi in plotindices
                dm = if @isdefined kinetics
                    kinetics[pi][i]
                else
                    damodels[i]
                end
                modelage(c, tT[pi], dm)
                plot!(h, c.tsteps_degassing, c.model_fraction, label="", color=:powderblue, alpha=0.15)
                plot!(h, c.tsteps_degassing[t], c.model_fraction[t], label="", color=:mediumblue, alpha=0.15)
            end
            plot!(c.tsteps_experimental, c.fraction_experimental, color=:black, label="", alpha=0.15, lw=2)
            plot!(c.tsteps_experimental[c.fit], c.fraction_experimental[c.fit], label="Data", color=:black, lw=2)
            plot!(Float64[], Float64[], label="Model", color=:mediumblue)
            savefig(h, "$(name)_$(c.name)_cumulative_degassing.pdf")
            display(h)
        end
    end

## --- Plot calculated and observed step heating age spectra

    for C in (SingleDomain, MultipleDomain)
        for i in eachindex(chrons)
            c = chrons[i]
            if c isa C  
                agedist = tT.resultdist[i,:]
                modelage = sort(agedist)
                modelfraction = range(0, 1, length(modelage))
                h = plot(modelfraction, modelage,
                    title = "$(c.name) $C step heating",
                    xlabel = "Fraction released",
                    ylabel = (eltype(c) <: Thermochron.HeliumSample) ? "Rstep/Rbulk" : "Age [Ma]",
                    label = "",
                    legend = :topleft,
                    framestyle = :box,
                    color = :powderblue,
                    lw = 2,
                )
                errorbox!(h, c, color=:black)
                fmin = (findfirst(c.fit) > 1) ? c.fraction_experimental[findfirst(c.fit)-1] : 0
                fmax = c.fraction_experimental[findlast(c.fit)]
                t = fmin .<= modelfraction .<= fmax
                plot!(h, modelfraction[t], modelage[t],
                    label = "Model (average)",
                    seriestype = :line,
                    color = :mediumblue,
                    lw = 2,
                )
                savefig(h, "$(name)_$(c.name)_step_heating.pdf")
                display(h)
            end
        end
    end

    if @isdefined kinetics
        # SDDiffusivity
        im = findall(x->isa(x, SDDiffusivity), damodels[:,1])
        id = findall(x->isa(x, SDDiffusivity), kinetics.dmdist[:,1])
        for i in eachindex(im,id)
            grain_name = chrons[im[i]].name
            hdm = plot(damodels[im[i]], kinetics.dmdist[id[i],:]; title = "$grain_name SDDiffusivity")
            savefig(hdm, "$(name)_$(grain_name)_SDDiffusivity_kinetics.pdf")
            display(hdm)
        end
        # MDDiffusivity
        im = findall(x->isa(x, MDDiffusivity), damodels[:,1])
        id = findall(x->isa(x, MDDiffusivity), kinetics.dmdist[:,1])
        for i in eachindex(im,id)
            r = last(first(chrons[im[i]].domains).redges)
            grain_name = chrons[im[i]].name
            hdm = plot(damodels[im[i]], kinetics.dmdist[id[i],:], r; title = "$grain_name MDDiffusivity")
            savefig(hdm, "$(name)_$(grain_name)_MDDiffusivity_kinetics.pdf")
            display(hdm)
        end
    end

## --- Create image of paths in the posterior distribution

    @time (tTimage, xc, yc) = image_from_paths!(tT; xresolution=1800, yresolution=1200, method=:nearest, xrange=boundary.agepoints, yrange=boundary.T₀)

## --- Plot image with 'ylcn' custom colorscale
 
    # Prepare axes
    k = plot(layout = grid(1,2, widths=[0.94, 0.06]), framestyle=:box)

    # Plot image with colorscale in first subplot
    colormap = ylcn # Type `colormaps` to see a list of available options
    A = imsc(tTimage, colormap, 0, nanpctile(tTimage[:],98.5))
    plot!(k[1], xlabel="Time (Ma)", ylabel="Temperature (°C)", legend=false, framestyle=:box, yflip=true, xflip=true, minorticks=true, tick_dir=:out,)
    plot!(k[1], xc, yc, A, aspectratio=params.tinit/params.Tinit/1.5, xlims=(0,params.tinit), ylims=(params.Tnow,params.Tinit))
    plot!(k[1], constraint, color=:red, lw=1) # Add constraint boxes, if any

    # Add colorbar in second subplot
    cb = imsc(repeat(0:100, 1, 10), colormap, 0, 100)
    plot!(k[2], 0:0.01:0.1, 0:0.01:1, cb, ylims=(0,1), xticks=false, tick_dir=:out, framestyle=:box, yflip=false, ylabel="Relative path density", guide_position=:right)

    savefig(k, name*"_tT.pdf")
    display(k)

## --- (Optional) Plot a zoomed-in version of the posterior path image

    # xrange = (0, 10)    # time [Ma]
    # yrange = (0, 400)   # Temperature [C]
    # @time (tTimageZoom, xcZoom, ycZoom) = image_from_paths!(tT; xresolution=1800, yresolution=1200, method=:nearest, xrange, yrange)

    # # Prepare axes
    # k = plot(layout = grid(1,2, widths=[0.94, 0.06]), framestyle=:box)

    # # Plot image with colorscale in first subplot
    # colormap = ylcn # Type `colormaps` to see a list of available options
    # A = imsc(tTimageZoom, colormap, 0, nanpctile(tTimage[:],98.5))
    # plot!(k[1], xlabel="Time (Ma)", ylabel="Temperature (°C)", legend=false, framestyle=:box, yflip=true, xflip=true, minorticks=true, tick_dir=:out,)
    # plot!(k[1], xcZoom, ycZoom, A, aspectratio=nanrange(xrange)/nanrange(yrange)/1.5, xlims=xrange, ylims=yrange)
    # plot!(k[1], constraint, color=:red, lw=1,) # Add constraint boxes, if any

    # # Add colorbar in second subplot
    # cb = imsc(repeat(0:100, 1, 10), colormap, 0, 100)
    # plot!(k[2], 0:0.01:0.1, 0:0.01:1, cb, ylims=(0,1), xticks=false, tick_dir=:out, framestyle=:box, yflip=false,  ylabel="Relative path density", guide_position=:right)

    # savefig(k, name*"_tT_$(xrange[1])_$(xrange[2])Ma_$(yrange[1])_$(yrange[2])C.pdf")
    # display(k)

## --- End of File

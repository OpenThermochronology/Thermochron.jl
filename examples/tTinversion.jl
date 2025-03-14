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
    # (12 ZirconHe, 5 ApatiteHe, 47 ApatiteFT, 269 ApatiteTrackLength)
    name = "Manitoba"
    ds = importdataset("manitoba.csv", ',', importas=:Tuple)

    # # Literature samples from Guenthner et al. 2013 (AJS), Minnesota
    # # (23 ZirconHe, 11 ApatiteHe)
    # name = "Minnesota"
    # ds = importdataset("minnesota.csv", ',', importas=:Tuple)

    # # OL13 multiple domain diffusion example
    # name = "OL13"
    # ds = importdataset("ol13.csv", ',', importas=:Tuple)

## --- Prepare problem

    dt = 8.0 # [Ma] Model timestep
    tinit = ceil(maximum(ds.crystallization_age_Ma)/dt) * dt # [Ma] Model start time

    model = (
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
        tsteps = cntr(0:dt:tinit),      # [Ma] Forward time discretiziation (from start of model)
        agesteps = cntr(tinit:-dt:0),   # [Ma] Age discretization (relative to the present)
        minpoints = 15,                 # [n] Minimum allowed number of model t-T points (nodes)
        maxpoints = 50,                 # [n] Maximum allowed number of model t-T points (nodes)
        dynamicsigma = false,           # Update model uncertainties throughout inversion?
        dynamicjumping = true,          # Update the t and T jumping (proposal) distributions based on previously accepted jumps
        # Damage and annealing models for diffusivity (specify custom kinetics if desired)
        adm = RDAAM(),                  # Flowers et al. 2009 (doi: 10.1016/j.gca.2009.01.015) apatite diffusivity model
        zdm = ZRDAAM(),                 # Guenthner et al. 2013 (doi: 10.2475/03.2013.01) zircon diffusivity model
        aftm = Ketcham1999FC(),         # Ketcham et al. 2007 (doi: 10.2138/am.2007.2281) apatite fission track model
        zftm = Yamada2007PC(),          # Yamada et al. 2007 (doi: 10.1016/j.chemgeo.2006.09.002) zircon fission track model
        mftm = Jones2021FA(),           # Re-fit from Jones et al. 2021 (doi: 10.5194/gchron-3-89-2021) monazite fission track 
        # Optional simulated annealing during burnin, wherein p_accept = max(exp((llₚ-ll)/T), 1)
        # T = T0annealing * exp(-λannealing * n) + 1 at step number n of burnin
        T0annealing = 5,                # [unitless] initial annealing "temperature" (set to 0 for no simulated annealing).
    )

    # Default: no detail interval
    detail = DetailInterval()

    # Uncomment this section to require greater t-T node density in some time interval
    # (typically the youngest end of the total time interval, where you may expect the data more resolving power)
    detail = DetailInterval(
        agemin = 0.0, # Youngest end of detail interval
        agemax = 1000, # Oldest end of detail interval
        minpoints = 8, # Minimum number of points in detail interval
    )

    # Boundary conditions (e.g. 10C at present and 650 C at the time of zircon formation).
    boundary = Boundary(
        agepoints = [0.0, model.tinit],          # [Ma] Final and initial time
        T₀ = [model.Tnow, model.Tinit],          # [C] Final and initial temperature
        ΔT = [model.ΔTnow, model.ΔTinit],        # [C] Final and initial temperature range (positive or negative)
        tboundary = :reflecting, # Reflecting time boundary conditions
        Tboundary = :reflecting, # Reflecting temperature boundary conditions
    )

    # Default: No constraints are imposed
    constraint = Constraint()

    # # Uncomment this section if you wish to impose an unconformity or other constraint
    # # at any point in the record.
    # constraint = Constraint(
    #     agedist = [ Normal(510,20),],  # [Ma] Age distribution
    #     Tdist =   [  Uniform(0,50),],  # [C] Temperature distribution
    # )
    # name *= "_constrained"

## --- Process data into Chronometer objects

    chrons, damodels = chronometers(ds, model)

    # Model uncertainty is not well known (depending on annealing parameters,
    # decay constants, diffusion parameters, etc.), but is certainly non-zero.
    # In addition, observed thermochronometric ages often display excess dispersion beyond analytical
    # uncertainty. Here we specify a default minimum uncertainty representing the average expected 
    # misfit between model and data, which may optionally be resampled during inversion.
    model = (;model...,
        σcalc = fill(25., length(chrons)),     # [Ma] model uncertainty
    )

## --- Age uncertainty resampling: estimate expected misfit (σcalc) from excess dispersion of the data itself

    # Empirical age uncertainty for apatite
    tap = isa.(chrons, ApatiteHe)
    if count(tap) > 1
        h = ageeuplot(chrons[tap], label="Internal uncertainty", title="apatite")
        empiricaluncertainty!(model.σcalc, chrons, ApatiteHe)
        σtotal = sqrt.(get_age_sigma(chrons[tap]).^2 + model.σcalc[tap].^2)
        ageeuplot!(h, chrons[tap], yerror=2*σtotal, label="Empirical uncertainty")
        display(h)
    end

    # Empirical age uncertainty for zircon
    tzr = isa.(chrons, ZirconHe)
    if count(tzr) > 1
        h = ageeuplot(chrons[tzr], label="Internal uncertainty", title="zircon")
        empiricaluncertainty!(model.σcalc, chrons, ZirconHe)
        σtotal = sqrt.(get_age_sigma(chrons[tzr]).^2 + model.σcalc[tzr].^2)
        ageeuplot!(h, chrons[tzr], yerror=2*σtotal, label="Empirical uncertainty")
        display(h)
    end

## --- Invert for maximum likelihood t-T path

    # Run Markov Chain
    # @time tT = MCMC(chrons, damodels, model, boundary, constraint, detail)
    @time tT, kinetics = MCMC_varkinetics(chrons, damodels, model, boundary, constraint, detail)
    @info """tT.tpointdist & tT.Tpointdist collected, size: $(size(tT.Tpointdist))
    Mean log-likelihood: $(nanmean(view(tT.lldist, model.burnin:model.nsteps)))
    Mean acceptance rate: $(nanmean(view(tT.acceptancedist, model.burnin:model.nsteps)))
    Mean npoints: $(nanmean(view(tT.ndist, model.burnin:model.nsteps)))
    Mean jₜ: $(nanmean(view(tT.jtdist,model.burnin:model.nsteps)))
    Mean jT: $(nanmean(view(tT.jTdist, model.burnin:model.nsteps)))
    """

    # # Save tTs using JLD (compressed)
    # using JLD: jldopen, @write
    # jldopen("$name.jld", "w", compress=true) do file
    #     @write file tT
    #     @write file model
    # end

    # # Alternatively, save as MAT file
    # using MAT
    # matwrite("$name.mat", Dict(
    #     "tpointdist"=>tT.tpointdist,
    #     "Tpointdist"=>tT.Tpointdist,
    #     "ndist"=>tT.ndist,
    #     "resultdist"=>tT.resultdist,
    #     "lldist"=>tT.lldist,
    #     "acceptancedist"=>tT.acceptancedist,
    #     "model"=>Dict(
    #         replace.(string.(keys(model)), "σ"=>"sigma", "λ"=>"lambda", "Δ"=>"Delta") .=> values(model)
    #     )
    # ), compress=true)

    # Plot log likelihood distribution
    h = plot(tT.lldist, xlabel="Step number", ylabel="Log likelihood", label="", framestyle=:box)
    savefig(h, name*"_lldist.pdf")
    display(h)

## --- Plot distribution of number of model t-T points (nodes)

    h = plot(tT.ndist, label="", framestyle=:box)
    plot!(xlabel="Step number", ylabel="Number of model t-T nodes", ylims=(0,last(ylims())))
    savefig(h, name*"_ndist.pdf")
    display(h)

## --- Plot moving average of acceptance distribution

    h = plot(movmean(tT.acceptancedist,100), label="", framestyle=:box)
    plot!(xlabel="Step number", ylabel="Acceptance probability (mean of 100)", ylims=(0,1))
    savefig(h, name*"_acceptance.pdf")
    display(h)

## --- Plot calculated/observed ages as a function of eU (ZirconHe, ApatiteHe)

    C = (ApatiteHe, ZirconHe,)
    mincolor = ("apatite", "zircon",)
    for i in eachindex(C, mincolor)
        t = isa.(chrons, C[i])
        if any(t)
            σtotal = sqrt.(get_age_sigma(chrons[t]).^2 + model.σcalc[t].^2)
            h = ageeuplot(chrons[t], yerror=2σtotal,
                label="Data (2σ total)", 
                color = :black, 
                title = "$(C[i])",
            )
            agedist = tT.resultdist[t,:]
            m = nanmean(agedist, dim=2)
            l = nanpctile(agedist, 2.5, dim=2)
            u = nanpctile(agedist, 97.5, dim=2)
            scatter!(h, eU.(chrons[t]), m, 
                yerror=(m-l, u-m), 
                label="Model (95%CI)", 
                color=mineralcolors[mincolor[i]], 
                msc=mineralcolors[mincolor[i]], 
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
                savefig(hdm, "$(name)_$(D)_kinetics.pdf")
                display(hdm)
            end
        end
    end

## -- Plot calculated/observed ages as a function of rmr0 (ApatiteFT)

    t = isa.(chrons, ApatiteFT)
    if any(t)
        rmr0 = chrons[t] .|> x->x.rmr0
        μobs = get_age(chrons[t])
        σtotal = sqrt.(get_age_sigma(chrons[t]).^2 + model.σcalc[t].^2)
        h = scatter(rmr0, μobs, yerror=2σtotal,
            xlabel = "rmr0 [unitless]",
            ylabel = "Age [Ma]",
            label = "Data (2σ total)", 
            framestyle = :box,
            color = :black,
            title = "ApatiteFT",
        )
        agedist = tT.resultdist[t,:]
        m = nanmean(agedist, dim=2)
        l = nanpctile(agedist, 2.5, dim=2)
        u = nanpctile(agedist, 97.5, dim=2)
        scatter!(h, rmr0, m, 
            yerror = (m-l, u-m), 
            label = "Model (95%CI)", 
            color = mineralcolors["apatite"], 
            msc = mineralcolors["apatite"], 
        )
        savefig(h, "$(name)_ApatiteFT_predicted.pdf")
        display(h)
    end

## -- Plot calculated and observed ages (most other chronometers)

    C = (PlanarAr, PlanarHe, SphericalAr, SphericalHe, ZirconFT, MonaziteFT)
    mincolor = ("feldspar", "hematite", "feldspar", "hematite", "zircon", "monazite")
    for i in eachindex(C, mincolor)
        t = isa.(chrons, C[i])
        if any(t)
            μobs = get_age(chrons[t])
            σtotal = sqrt.(get_age_sigma(chrons[t]).^2 + model.σcalc[t].^2)
            h = plot(μobs, yerror=2σtotal,
                xlabel = "Sample number",
                ylabel = "Age [Ma]",
                label = "Data (2σ total)", 
                framestyle = :box,
                color = :black,
                title = "$(C[i])",
            )
            agedist = tT.resultdist[t,:]
            m = nanmean(agedist, dim=2)
            l = nanpctile(agedist, 2.5, dim=2)
            u = nanpctile(agedist, 97.5, dim=2)
            scatter!(h, m, 
                yerror = (m-l, u-m), 
                label = "Model (95%CI)", 
                color = mineralcolors[mincolor[i]], 
                msc = mineralcolors[mincolor[i]], 
            )
            savefig(h, "$(name)_$(C[i])_predicted.pdf")
            display(h)
        end
    end

## -- Fission track length histograms (ZirconTrackLength, ApatiteTrackLength, MonaziteTrackLength)

    C = (ZirconTrackLength, ApatiteTrackLength, MonaziteTrackLength)
    mincolor = ("zircon", "apatite", "monazite")
    for i in eachindex(C, mincolor)
        t = isa.(chrons, C[i])
        if any(t)
            h = histogram(Thermochron.val.(chrons[t]), bins=0:0.25:20, 
                normalized=true,
                xlims = (0,20),
                xlabel = "Track length [μm]",
                ylabel = "Probability density",
                label = "Data (N=$(count(t)))", 
                framestyle = :box,
                legend = :topleft,
                color = :black,
                alpha = 0.75,
                title = "$(C[i])",
            )
            lengthdist = tT.resultdist[t,:]
            histogram!(h, vec(lengthdist), bins=0:0.25:20, 
                normalized=true, 
                label = "Model",
                color = mineralcolors[mincolor[i]],
                fill = true,
                alpha = 0.75,
            )
            savefig(h, "$(name)_$(C[i])_predicted.pdf")
            display(h)
        end
    end

## --- Plot calculated and observed MDD age spectra

    for i in eachindex(chrons)
        c = chrons[i]
        if c isa MultipleDomain  
            agedist = tT.resultdist[i,:]
            modelage = sort(agedist)
            modelfraction = range(0, 1, length(modelage))
            h = plot(modelfraction, modelage,
                seriestype = :line,
                color = :powderblue,
                framestyle = :box,
                legend = :topleft,
                xlabel = "Fraction released",
                ylabel = "Age [Ma]",
                label = "",
            )
            errorbox!(h, c, color=:black)
            t = minimum(c.fraction_experimental[c.fit]) .<= modelfraction .<= maximum(c.fraction_experimental[c.fit])
            plot!(h, modelfraction[t], modelage[t],
                seriestype=:line,
                color = :mediumblue,
                lw = 2,
                label = "Model (average)",
            )
            savefig(h, "$(name)_$(ds.grain_name[i])_predicted.pdf")
            display(h)
        end
    end

    if @isdefined kinetics
        im = findall(x->isa(x, MDDiffusivity), damodels[:,1])
        id = findall(x->isa(x, MDDiffusivity), kinetics.dmdist[:,1])
        for i in eachindex(im,id)
            r = last(first(chrons[im[i]].domains).redges)
            hdm = plot(damodels[im[i]], kinetics.dmdist[id[i],:], r)
            savefig(hdm, "$(name)_$(ds.grain_name[im[i]])_MDD_kinetics.pdf")
            display(hdm)
        end
    end

## --- Create image of paths

    @time (tTimage, xc, yc) = image_from_paths!(tT; xresolution=1800, yresolution=1200, xrange=boundary.agepoints, yrange=boundary.T₀)

## --- Plot image with 'ylcn' custom colorscale
 
    # Prepare axes
    k = plot(layout = grid(1,2, widths=[0.94, 0.06]), framestyle=:box)

    # Plot image with colorscale in first subplot
    A = imsc(tTimage, ylcn, 0, nanpctile(tTimage[:],98.5))
    plot!(k[1], xlabel="Time (Ma)", ylabel="Temperature (°C)", yticks=model.Tnow:50:model.Tinit, xticks=model.tnow:500:model.tinit, yminorticks=5, xminorticks=5, tick_dir=:out, framestyle=:box)
    plot!(k[1], xc, yc, A, yflip=true, xflip=true, legend=false, aspectratio=model.tinit/model.Tinit/1.5, xlims=(0,model.tinit), ylims=(model.Tnow,model.Tinit))
    plot!(k[1], constraint) # Add constraint boxes

    # Add colorbar in second subplot
    cb = imsc(repeat(0:100, 1, 10), ylcn, 0, 100)
    plot!(k[2], 0:0.01:0.1, 0:0.01:1, cb, ylims=(0,1), xticks=false, framestyle=:box, yflip=false)

    savefig(k, name*"_tT.pdf")
    display(k)

## --- Plot a zoomed-in version

    # xrange = (0,100)
    # yrange=(-10,60)
    # @time (tTimageZoom, xcZoom, ycZoom) = image_from_paths!(tT; xresolution=1800, yresolution=1200, xrange, yrange)

    # # Prepare axes
    # k = plot(layout = grid(1,2, widths=[0.94, 0.06]), framestyle=:box)

    # # Plot image with colorscale in first subplot
    # A = imsc(tTimageZoom, ylcn, 0, nanpctile(tTimage[:],98.5))
    # plot!(k[1], xlabel="Time (Ma)", ylabel="Temperature (°C)", tick_dir=:out, framestyle=:box)
    # plot!(k[1], xcZoom, ycZoom, A, yflip=true, xflip=true, legend=false, xlims=xrange, ylims=yrange)

    # # Add colorbar in second subplot
    # cb = imsc(repeat(0:100, 1, 10), ylcn, 0, 100)
    # plot!(k[2], 0:0.01:0.1, 0:0.01:1, cb, ylims=(0,1), xticks=false, framestyle=:box, yflip=false)

    # savefig(k, name*"_tT.pdf")
    # display(k)

## --- End of File
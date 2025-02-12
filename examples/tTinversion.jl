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
#   © 2024 C. Brenhin Keller and Kalin McDannell                                #                                                                            #
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

    # # Literature samples from Guenthner et al. 2013 (AJS), Minnesota
    # # (23 ZirconHe, 11 ApatiteHe)
    # name = "Minnesota"
    # ds = importdataset("minnesota.csv", ',', importas=:Tuple)

    # Literature samples from McDannell et al. 2022 (doi: 10.1130/G50315.1), Manitoba
    # (12 ZirconHe, 5 ApatiteHe, 47 ApatiteFT, 269 ApatiteTrackLength)
    name = "Manitoba"
    ds = importdataset("manitoba.csv", ',', importas=:Tuple)


## --- Prepare problem

    dt = 8.0 # [Ma] Model timestep
    tinit = ceil(maximum(ds.crystallization_age_Ma)/dt) * dt # [Ma] Model start time

    model = (
        nsteps = 40000,                 # [n] How many steps of the Markov chain should we run?
        burnin = 10000,                 # [n] How long should we wait for MC to converge (become stationary)
        dr = 1.0,                       # [μm] Radius step size
        dTmax = 10.0,                   # [Ma/dt] Maximum reheating/burial per model timestep. If too high, may cause numerical problems in Crank-Nicholson solve
        Tinit = 400.0,                  # [C] Initial model temperature (i.e., crystallization temperature)
        ΔTinit = -100.0,                # [C] Tinit can vary from Tinit to Tinit+ΔTinit
        Tnow = 0.0,                     # [C] Current surface temperature
        ΔTnow = 20.0,                   # [C] Tnow may vary from Tnow to Tnow+ΔTnow
        dt = dt,                        # [Ma] Model timestep
        tinit = tinit,                  # [Ma] Model start time
        tnow = 0.0,                     # [Ma] Model end time (today)
        tsteps = cntr(0:dt:tinit),      # [Ma] Forward time discretiziation (from start of model)
        agesteps = cntr(tinit:-dt:0),   # [Ma] Age discretization (relative to the present)
        minpoints = 15,                 # [n] Minimum allowed number of t-T points
        maxpoints = 50,                 # [n] Maximum allowed number of t-T points
        rescale = false,                # Attempt to hedge against systematic errors by limiting the log likeihood contribution of each chronometer to scale as sqrt(n) instead of n
        trackhist = false,
        dynamicsigma = false,           # Update model uncertainties throughout inversion?
        dynamicjumping = true,          # Update the t and t jumping (proposal) distributions based on previously accepted jumps
        # Damage and annealing models for diffusivity (specify custom kinetics if desired)
        adm = RDAAM(),                  # Flowers et al. 2009 (doi: 10.1016/j.gca.2009.01.015) apatite diffusivity model
        zdm = ZRDAAM(),                 # Guenthner et al. 2013 (doi: 10.2475/03.2013.01) zircon diffusivity model
        aftm = Ketcham1999FC(),         # Ketcham et al. 2007 (doi: 10.2138/am.2007.2281) apatite fission track model
        zftm = Yamada2007PC(),          # Yamada et al. 2007 (doi: 10.1016/j.chemgeo.2006.09.002) zircon fission track model
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
        agemax = 541.0, # Oldest end of detail interval
        minpoints = 7, # Minimum number of points in detail interval
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
    #     agedist = [ Normal(974,122),  Uniform(497,509),   Uniform(305,310)],  # [Ma] Age distribution
    #     Tdist =   [  Uniform(0,200),     Uniform(0,50),      Uniform(0,50)],  # [C] Temperature distribution
    # )
    # name *= "_constrained"

## --- Process data into Chronometer objects

    chrons = chronometers(ds, model)

    # Model uncertainty is not well known (depends on annealing parameters,
    # decay constants, diffusion parameters, etc.), but is certainly non-zero.
    # In addition, observed thermochronometric ages often display excess dispersion beyond analytical
    # uncertainty. Here we specify a default minimum uncertainty representing the average expected 
    # misfit between model and data, which may optionally be resampled during inversion.
    model = (;model...,
        σcalc = fill(25., length(chrons)),     # [Ma] model uncertainty
    )

## --- Age uncertainty resampling: estimate expected misfit (σcalc) from excess dispersion of the data itself

    # # Empirical age uncertainty for apatite
    # tap = isa.(chrons, ApatiteHe)
    # h = ageeuplot(chrons[tap], label="Internal uncertainty", title="apatite")
    # empiricaluncertainty!(model.σcalc, chrons, ApatiteHe)
    # σtotal = sqrt.(get_age_sigma(chrons[tap]).^2 + model.σcalc[tap].^2)
    # ageeuplot!(h, chrons[tap], yerror=2*σtotal, label="Empirical uncertainty")
    # display(h)

    # # Empirical age uncertainty for zircon
    # tzr = isa.(chrons, ZirconHe)
    # h = ageeuplot(chrons[tzr], label="Internal uncertainty", title="zircon")
    # empiricaluncertainty!(model.σcalc, chrons, ZirconHe)
    # σtotal = sqrt.(get_age_sigma(chrons[tzr]).^2 + model.σcalc[tzr].^2)
    # ageeuplot!(h, chrons[tzr], yerror=2*σtotal, label="Empirical uncertainty")
    # display(h)

## --- Invert for maximum likelihood t-T path

    # Run Markov Chain
    # @time tT = MCMC(chrons, model, boundary, constraint, detail)
    @time tT, kinetics = MCMC_varkinetics(chrons, model, boundary, constraint, detail)
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
    savefig(h, name*"_tT.lldist.pdf")
    display(h)

## --- Plot model ages vs observed ages in age-eU space (zircon)
    tz = isa.(chrons, ZirconHe)
    if any(tz)
        h = ageeuplot(chrons[tz],
            label="Data (2σ)", 
            color=:black, 
        )
        zircon_agedist = tT.resultdist[tz, :]
        m = nanmean(zircon_agedist, dims=2)
        l = nanpctile(zircon_agedist, 2.5, dims=2)
        u = nanpctile(zircon_agedist, 97.5, dims=2)
        scatter!(h, eU.(chrons[tz]), m, 
            yerror=(m-l, u-m), 
            label="Model + 95%CI", 
            color=mineralcolors["zircon"], 
            msc=mineralcolors["zircon"],
        )
        savefig(h, name*"_zircon_Age-eU.pdf")
        display(h)
    end

## --- Plot model ages vs observed ages in age-eU space (apatite)

    ta = isa.(chrons, ApatiteHe)
    if any(ta)
        h = ageeuplot(chrons[ta], 
            label="Data (2σ)", 
            color=:black, 
        )
        apatite_agedist = tT.resultdist[ta,:]
        m = nanmean(apatite_agedist, dims=2)
        l = nanpctile(apatite_agedist, 2.5, dims=2)
        u = nanpctile(apatite_agedist, 97.5, dims=2)
        scatter!(h, eU.(chrons[ta]), m, 
            yerror=(m-l, u-m), 
            label="Model + 95%CI", 
            color=mineralcolors["apatite"], 
            msc=mineralcolors["apatite"], 
        )
        savefig(h, name*"_apatite_Age-eU.pdf")
        display(h)
    end

## --- Plot moving average of acceptance distribution

    h = plot(movmean(tT.acceptancedist,100), label="", framestyle=:box)
    plot!(xlabel="Step number", ylabel="acceptance probability (mean of 100)")
    savefig(h, name*"_acceptance.pdf")
    display(h)

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

## --- End of File

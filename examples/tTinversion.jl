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

    # Literature samples from Guenthner et al. 2013 (AJS), Minnesota
    name = "Minnesota"
    ds = importdataset("minnesota.csv", ',', importas=:Tuple)

## --- Prepare problem

    model = (
        nsteps = 600000,            # How many steps of the Markov chain should we run?
        burnin = 350000,            # How long should we wait for MC to converge (become stationary)
        dr = 1.0,                   # Radius step, in microns
        dt = 8.0,                   # Time step size in Myr
        dTmax = 10.0,               # Maximum reheating/burial per model timestep. If too high, may cause numerical problems in Crank-Nicholson solve
        Tinit = 300.0,              # Initial model temperature (in C) (i.e., crystallization temperature)
        ΔTinit = -50.0,             # Tinit can vary from Tinit to Tinit+ΔTinit
        Tnow = 0.0,                 # Current surface temperature (in C)
        ΔTnow = 20.0,               # Tnow may vary from Tnow to Tnow+ΔTnow
        tnow = 0.0,                 # Today
        minpoints = 15,             # Minimum allowed number of t-T points
        maxpoints = 50,             # Maximum allowed number of t-T points
        dynamicsigma = true,        # Update model uncertainties?
        dynamicjumping = true,      # Update the t and t jumping (proposal) distributions based on previously accepted jumps
        # Damage and annealing models for diffusivity (specify custom kinetics if desired)
        adm = RDAAM(),
        zdm = ZRDAAM(), 
        aftm = SimplifiedCurvilinear(),
        # Model uncertainty is not well known (depends on annealing parameters,
        # decay constants, diffusion parameters, etc.), but is certainly non-zero.
        # Here we add (in quadrature) a blanket model uncertainty of 5 Ma.
        σmodel = 5.0,               # [Ma]
        σannealing = 35.0,          # initial annealing uncertainty [Ma]
        λannealing = 1 ./ 100_000   # annealing decay [1/n]
    )

    # Crystallization ages and start time
    tinit = ceil(maximum(ds.crystallization_age_Ma)/model.dt) * model.dt
    model = (model...,
        tinit = tinit,
        agesteps = Array{Float64}(tinit-model.dt/2 : -model.dt : 0+model.dt/2),
        tsteps = Array{Float64}(0+model.dt/2 : model.dt : tinit-model.dt/2),
    )

    # Default: no detail interval
    detail = DetailInterval()

    # # Uncomment this section to require greater t-T node density in some time interval
    # detail = DetailInterval(
    #     agemin = 0.0, # Youngest end of detail interval
    #     agemax = 541.0, # Oldest end of detail interval
    #     minpoints = 7, # Minimum number of points in detail interval
    # )

    # Boundary conditions (e.g. 10C at present and 650 C at the time of zircon formation).
    boundary = Boundary(
        agepoints = [model.tnow, model.tinit],   # [Ma] Final and initial time
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
    #     agedist = [Uniform(500,580),],  # [Ma] Age distribution
    #     Tdist =   [   Uniform(0,50),],  # [C] Temperature distribution
    # )
    # name *= "_unconf"

## --- Process data into Chronometer objects

    chrons = chronometers(ds, model)

## --- Age uncertainty resampling

    # Empirical age uncertainty for apatite
    tap = isa.(chrons, ApatiteHe)
    h = ageeuplot(chrons[tap], label="Internal uncertainty", title="apatite")
    empiricaluncertainty!(chrons, ApatiteHe)
    ageeuplot!(h, chrons[tap], label="Empirical uncertainty")
    display(h)

    # Empirical age uncertainty for zircon
    tzr = isa.(chrons, ZirconHe)
    h = ageeuplot(chrons[tzr], label="Internal uncertainty", title="zircon")
    empiricaluncertainty!(chrons, ZirconHe)
    ageeuplot!(h, chrons[tzr], label="Empirical uncertainty")
    display(h)


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

    # Desired rsolution of tTing image
    xresolution = 2000
    yresolution = 1000
    burnin = model.burnin

    # Resize the post-burnin part of the stationary distribution
    tTdist = Array{Float64}(undef, xresolution, model.nsteps-burnin)
    xq = range(0, model.tinit, length=xresolution)
    @time @inbounds for i = 1:model.nsteps-burnin
        linterp1s!(view(tTdist,:,i), view(tT.tpointdist,:,i+burnin), view(tT.Tpointdist,:,i+burnin), xq)
    end

    # Calculate composite image
    ybinedges = range(model.Tnow, model.Tinit, length=yresolution+1)
    tTimage = zeros(yresolution, size(tTdist,1))
    @time @inbounds for i=1:size(tTdist,1)
        histcounts!(view(tTimage,:,i), view(tTdist,i,:), ybinedges)
    end

## --- Plot image with 'ylcn' custom colorscale

    # Prepare axes
    k = plot(layout = grid(1,2, widths=[0.94, 0.06]), framestyle=:box)

    # Plot image with colorscale in first subplot
    A = imsc(tTimage, ylcn, 0, nanpctile(tTimage[:],98.5))
    plot!(k[1], xlabel="Time (Ma)", ylabel="Temperature (°C)", yticks=model.Tnow:50:model.Tinit, xticks=model.tnow:500:model.tinit, yminorticks=5, xminorticks=5, tick_dir=:out, framestyle=:box)
    plot!(k[1], xq, cntr(ybinedges), A, yflip=true, xflip=true, legend=false, aspectratio=model.tinit/model.Tinit/1.5, xlims=(0,model.tinit), ylims=(model.Tnow,model.Tinit))
    plot!(k[1], constraint) # Add constraint boxes

    # Add colorbar in second subplot
    cb = imsc(repeat(0:100, 1, 10), ylcn, 0, 100)
    plot!(k[2], 0:0.01:0.1, 0:0.01:1, cb, ylims=(0,1), xticks=false, framestyle=:box, yflip=false)

    savefig(k, name*"_tT.pdf")
    display(k)

## --- End of File

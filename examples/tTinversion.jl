# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                              tTinversion.jl                                   #
#                                                                               #
#       A transdimensional Bayesian Markov-chain Monte Carlo code for           #
#   inverting time-temperature paths from zircon U-Th/He ages, using a 1D       #
#   Crank-Nicholson finite difference diffusion model along with the zircon     #
#   damage and annealing model "ZRDAAM" from Guenthner et al, 2013 (AJS).       #
#   Requires uncorrected ZHe ages, grain size, [U] and [Th] in ppm.             #
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

    # Populate data NamedTuple from imported dataset
    data = (;
        halfwidth = copy(ds.halfwidth_um),            # Crystal half-width, in microns
        U = copy(ds.U238_ppm),                        # U concentration, in PPM
        Th = copy(ds.Th232_ppm),                      # Th-232 concentration, in PPM
        Sm = copy(ds.Sm147_ppm),                      # Sm-147 concentration, in PPM (optional)
        HeAge = copy(ds.raw_He_age_Ma),               # He age, in Ma
        HeAge_sigma = copy(ds.raw_He_age_sigma_Ma),   # He age uncertainty (1-sigma), in Ma
        crystAge = copy(ds.crystallization_age_Ma),   # Crystallization age, in Ma
        mineral = copy(ds.mineral),                   # i.e., "zircon" or "apatite"
    )

    ta = containsi.(data.mineral, "apatite")
    tz = containsi.(data.mineral, "zircon")
    eU = data.U+0.238*data.Th; # Used for plotting, etc.

## --- Empirical uncertainty estimation (using curve fit)

    # using LsqFit: curve_fit
    # function bilinearexponential(x,p)
    #     @. p[3:5] = abs(p[3:5])
    #     A, loc, scl, shp, skw = p
    #     xs = @. (x - loc)/scl # X standardized by location (mean-like) and scale (variance-like)
    #     v = @. 1/2 - atan(xs)/3.141592653589793 # Sigmoid (positive on LHS)
    #     @. exp(A + shp*skw*xs*v - shp/skw*xs*(1-v))
    # end
    # age_sigma = copy(ds.HeAge_Ma_sigma_raw)
    # age_sigma_empirical = similar(age_sigma)

    # if any(ta)
    #     p = [1, nanmean(eU[ta]), nanstd(eU[ta]), 1, 1]
    #     fobj = curve_fit(bilinearexponential,eU[ta],data.HeAge[ta],p)
    #     σ_external = nanstd(data.HeAge[ta] .- bilinearexponential(eU[ta],fobj.param))
    #     σ_external /= sqrt(2) # Assume half of variance is unknown external uncertainty
    #     @. age_sigma_empirical[ta] = sqrt(data.HeAge_sigma[ta]^2 + σ_external^2)

    #     h = plot(xlabel="eU", ylabel="Age", framestyle=:box, title="apatite")
    #     x = range(extrema(eU[ta])..., length=100)
    #     plot!(x, bilinearexponential(x,fobj.param), label="trendline")
    #     plot!(eU[ta], data.HeAge[ta], yerror=age_sigma_empirical[ta], seriestype=:scatter, c=:black, msc=:black, label="empirical")
    #     plot!(eU[ta], data.HeAge[ta], yerror=age_sigma[ta], seriestype=:scatter, c=mineralcolors["zircon"], msc=mineralcolors["zircon"], label="internal")
    #     display(h)
    # end

    # if any(tz)
    #     p = [1, nanmean(eU[tz]), nanstd(eU[tz]), 1, 1]
    #     fobj = curve_fit(bilinearexponential,eU[tz],data.HeAge[tz],p)
    #     σ_external = nanstd(data.HeAge[tz] .- bilinearexponential(eU[tz],fobj.param))
    #     σ_external /= sqrt(2) # Assume half of variance is unknown external uncertainty
    #     @. age_sigma_empirical[tz] = sqrt(data.HeAge_sigma[tz]^2 + σ_external^2)

    #     h = plot(xlabel="eU", ylabel="Age", framestyle=:box, title="zircon")
    #     x = range(extrema(eU[tz])..., length=100)
    #     plot!(x, bilinearexponential(x,fobj.param), label="trendline")
    #     plot!(eU[tz], data.HeAge[tz], yerror=age_sigma_empirical[tz], seriestype=:scatter, c=:black, msc=:black, label="empirical")
    #     plot!(eU[tz], data.HeAge[tz], yerror=age_sigma[tz], seriestype=:scatter, c=mineralcolors["zircon"], msc=mineralcolors["zircon"], label="internal")
    #     display(h)
    # end

## --- Empirical uncertainty estimation (using Gaussian window)

    age_sigma = copy(ds.HeAge_Ma_sigma_raw)
    age_sigma_empirical = similar(age_sigma)

    if any(ta)
        # Standard deviation of a Gaussian kernel in eU space, representing the 
        # range of eU over which zircons with similar eU should have similar ages
        σeU = 10

        # Calculate errors
        for i ∈ findall(ta)
            nearesteU = minimum(x->abs(x-eU[i]), eU[setdiff(findall(ta), i)])
            W = normpdf.(eU[i], max(σeU, nearesteU/2), eU[ta])
            σ_external = nanstd(data.HeAge[ta], W) # Weighted standard deviation
            σ_external /= sqrt(2) # Assume half of variance is unknown external uncertainty
            σ_internal = age_sigma[i]
            age_sigma_empirical[i] = sqrt(σ_external^2 + σ_internal^2)
        end

        h = plot(xlabel="eU", ylabel="Age", framestyle=:box, title="apatite")
        plot!(eU[ta], data.HeAge[ta], yerror=age_sigma_empirical[ta], seriestype=:scatter, c=:black, msc=:black, label="empirical")
        plot!(eU[ta], data.HeAge[ta], yerror=age_sigma[ta], seriestype=:scatter, c=mineralcolors["apatite"], msc=mineralcolors["apatite"], label="internal")
        display(h)
    end

    if any(tz)
        # Standard deviation of a Gaussian kernel in eU space, representing the 
        # range of eU over which zircons with similar eU should have similar ages
        σeU = 100

        # Calculate errors
        for i ∈ findall(tz)
            nearesteU = minimum(x->abs(x-eU[i]), eU[setdiff(findall(tz), i)])
            W = normpdf.(eU[i], max(σeU, nearesteU/2), eU[tz])
            σ_external = nanstd(data.HeAge[tz], W) # Weighted standard deviation
            σ_external /= sqrt(2) # Assume half of variance is unknown external uncertainty
            σ_internal = age_sigma[i]
            age_sigma_empirical[i] = sqrt(σ_external^2 + σ_internal^2)
        end

        h = plot(xlabel="eU", ylabel="Age", framestyle=:box, title="zircon")
        plot!(eU[tz], data.HeAge[tz], yerror=age_sigma_empirical[tz], seriestype=:scatter, c=:black, msc=:black, label="empirical")
        plot!(eU[tz], data.HeAge[tz], yerror=age_sigma[tz], seriestype=:scatter, c=mineralcolors["zircon"], msc=mineralcolors["zircon"], label="internal")
        display(h)
    end

## --- Prepare problem

    # Use empirical ages for zircon, 10 % for apatite
    data.HeAge_sigma[tz] .= age_sigma_empirical[tz]
    data.HeAge_sigma[ta] .= age_sigma_empirical[ta]

    model = (
        nsteps = 600_000,         # How many steps of the Markov chain should we run?
        burnin = 350_000,           # How long should we wait for MC to converge (become stationary)
        dr = 1.0,                   # Radius step, in microns
        dt = 10.0,                  # Time step size in Myr
        dTmax = 10.0,               # Maximum reheating/burial per model timestep. If too high, may cause numerical problems in Crank-Nicholson solve
        Tinit = 400.0,              # Initial model temperature (in C) (i.e., crystallization temperature)
        ΔTinit = -50.0,             # Tinit can vary from Tinit to Tinit+ΔTinit
        Tnow = 0.0,                 # Current surface temperature (in C)
        ΔTnow = 20.0,               # Tnow may vary from Tnow to Tnow+ΔTnow
        tnow = 0.0,                 # Today
        tinitMax = 3500.0,          # Ma -- forbid anything older than this
        minpoints = 20,             # Minimum allowed number of t-T points
        maxpoints = 50,             # Maximum allowed number of t-T points
        simplified = false,         # Prefer simpler tT paths?
        dynamicjumping = true,      # Update the t and t jumping (proposal) distributions based on previously accepted jumps
        # Damage and annealing models for diffusivity (specify custom kinetics if desired)
        adm = RDAAM(),
        zdm = ZRDAAM(), 
        # Model uncertainty is not well known (depends on annealing parameters,
        # decay constants, diffusion parameters, etc.), but is certainly non-zero.
        # Here we add (in quadrature) a blanket model uncertainty of 25 Ma.
        σmodel = 0.0,              # [Ma]
        σannealing = 135.0,          # initial annealing uncertainty [Ma]
        λannealing = 10 ./ 100_000  # annealing decay [1/n]
    )

    # Sort out crystallization ages and start time
    map!(x->min(x, model.tinitMax), data.crystAge, data.crystAge)
    tinit = ceil(maximum(data.crystAge)/model.dt) * model.dt
    model = (model...,
        tinit = tinit,
        agesteps = Array{Float64}(tinit-model.dt/2 : -model.dt : 0+model.dt/2),
        tsteps = Array{Float64}(0+model.dt/2 : model.dt : tinit-model.dt/2),
    )

    # Default: no detail interval
    detail = DetailInterval()

    # Uncomment this section to require greater t-T node density in some time interval
    detail = DetailInterval(
        agemin = 0.0, # Youngest end of detail interval
        agemax = 541.0, # Oldest end of detail interval
        minpoints = 7, # Minimum number of points in detail interval
    )

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

## --- Invert for maximum likelihood t-T path

    # Run Markov Chain
    # @time tT = MCMC(data, model, boundary, constraint, detail)
    @time tT, kinetics = MCMC_varkinetics(data, model, boundary, constraint, detail)
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
    #     "HeAgedist"=>tT.HeAgedist,
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

    if any(tz)
        h = scatter(eU[tz], data.HeAge[tz], 
            yerror=2*data.HeAge_sigma[tz], 
            label="Data (2σ)", 
            color=:black, 
            xlabel="eU (ppm)",
            ylabel="Age (Ma)",
            framestyle=:box,
        )
        zircon_agedist = tT.HeAgedist[tz, model.burnin:end]
        m = nanmean(zircon_agedist, dims=2)
        l = nanpctile(zircon_agedist, 2.5, dims=2)
        u = nanpctile(zircon_agedist, 97.5, dims=2)
        scatter!(h, eU[tz], m, 
            yerror=(m-l, u-m), 
            label="Model + 95%CI", 
            color=mineralcolors["zircon"], 
            msc=mineralcolors["zircon"],
        )
        savefig(h, name*"_zircon_Age-eU.pdf")
        display(h)
    end

## --- Plot model ages vs observed ages in age-eU space (apatite)

    if any(ta)
        h = scatter(eU[ta], data.HeAge[ta], 
            yerror=2*data.HeAge_sigma[ta], 
            label="Data (2σ)", 
            color=:black, 
            xlabel="eU (ppm)",
            ylabel="Age (Ma)",
            framestyle=:box,
        )
        apatite_agedist = tT.HeAgedist[ta,model.burnin:end]
        m = nanmean(apatite_agedist, dims=2)
        l = nanpctile(apatite_agedist, 2.5, dims=2)
        u = nanpctile(apatite_agedist, 97.5, dims=2)
        scatter!(h, eU[ta], m, 
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
    plot!(k[1], xlabel="Time (Ma)", ylabel="Temperature (°C)", yticks=0:50:400, xticks=0:500:3500, yminorticks=5, xminorticks=5, tick_dir=:out, framestyle=:box)
    plot!(k[1], xq, cntr(ybinedges), A, yflip=true, xflip=true, legend=false, aspectratio=model.tinit/model.Tinit/1.5, xlims=(0,model.tinit), ylims=(0,400))

    # Add colorbar in second subplot
    cb = imsc(repeat(0:100, 1, 10), ylcn, 0, 100)
    plot!(k[2], 0:0.01:0.1, 0:0.01:1, cb, ylims=(0,1), xticks=false, framestyle=:box, yflip=false)

    plot!([635.5, 717.4, 717.4, 635.5, 635.5],[0, 0, 650, 650, 0], fill=true, color=:white, alpha=0.5) #Sturtian & Marinoan glacial

    savefig(k, name*"_tT.pdf")
    display(k)

## ---

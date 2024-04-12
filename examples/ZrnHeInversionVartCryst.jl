# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                            ZrnHeInversion.jl                                  #
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
#   © 2022 C. Brenhin Keller and Kalin McDannell                                #                                                                            #
#                                                                               #
#       If running for the first time, you might consider instantiating the     #
#   manifest that came with this file                                           #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## ---  Load required packages

    using Thermochron
    using StatGeochem, Plots

    using LinearAlgebra
    # Diminishing returns with more than ~4 threads
    BLAS.get_num_threads() > 4 && BLAS.set_num_threads(4)

    # Make sure we're running in the directory where the script is located
    cd(@__DIR__)

    # # # # # # # # # # Choice of regional thermochron data # # # # # # # # # #

    # Literature samples from Guenthner et al. 2013 (AJS), Minnesota
    name = "Minnesota"
    ds = importdataset("minnesota.csv", ',', importas=:Tuple)

    # Populate data NamedTuple from imported dataset
    data = (
        halfwidth = ds.Halfwidth_um,            # Crystal half-width, in microns
        U = ds.U238_ppm,                        # U concentration, in PPM
        Th = ds.Th232_ppm,                      # Th-232 concentration, in PPM
        Sm = ds.Sm147_ppm,                      # Sm-147 concentration, in PPM (optional)
        HeAge = ds.HeAge_Ma_raw,                # He age, in Ma
        HeAge_sigma = ds.HeAge_Ma_sigma_10pct,  # He age uncertainty (1-sigma), in Ma
        crystAge = ds.CrystAge_Ma,              # Crystallization age, in Ma
        mineral = ds.Mineral
    );

## --- Prepare problem

    model = (
        nsteps = 20_000, # How many steps of the Markov chain should we run?
        burnin = 10_000, # How long should we wait for MC to converge (become stationary)
        dr = 1.0,    # Radius step, in microns
        dt = 10.0,   # Time step size in Myr
        dTmax = 25.0, # Maximum reheating/burial per model timestep. If too high, may cause numerical problems in Crank-Nicholson solve
        Tinit = 400.0, # Initial model temperature (in C) (i.e., crystallization temperature)
        ΔTinit = -50.0, # Tinit can vary from Tinit to Tinit+ΔTinit
        Tnow = 0.0, # Current surface temperature (in C)
        ΔTnow = 20.0, # Tnow may vary from Tnow to Tnow+ΔTnow
        tinitMax = 3500.0, # Ma -- forbid anything older than this
        minpoints = 20, # Minimum allowed number of t-T points
        maxpoints = 50, # Maximum allowed number of t-T points
        simplified = false, # Prefer simpler tT paths?
        dynamicjumping = true, # Update the t and t jumping (proposal) distributions based on previously accepted jumps
        # Diffusion parameters (ZRDAAM standard)
        DzEa = 165.0, # kJ/mol
        DzD0 = 193188.0, # cm^2/sec
        DN17Ea = 71.0, # kJ/mol
        DN17D0 = 6.367E-3, # cm^2/sec
        # Model uncertainty is not well known (depends on annealing parameters,
        # decay constants, diffusion parameters, etc.), but is certainly non-zero.
        # Here we add (in quadrature) a blanket model uncertainty of 25 Ma.
        σmodel = 25.0, # Ma
        σannealing = 35.0, # initial annealing uncertainty [Ma]
        λannealing = 10 ./ 100_000 # annealing decay [1/n]
    )

    # Sort out crystallization ages and start time
    map!(x->max(x, model.tinitMax), data.crystAge, data.crystAge)
    tinit = ceil(maximum(data.crystAge)/model.dt) * model.dt
    model = (model...,
        tinit = tinit,
        agesteps = Array{Float64}(tinit-model.dt/2 : -model.dt : 0+model.dt/2),
        tsteps = Array{Float64}(0+model.dt/2 : model.dt : tinit-model.dt/2),
    )

    detail = DetailInterval(
        agemin = 0.0, # Youngest end of detail interval
        agemax = 541.0, # Oldest end of detail interval
        minpoints = 7, # Minimum number of points in detail interval
    )

    # Boundary conditions (e.g. 10C at present and 650 C at the time of zircon formation).
    boundary = Boundary(
        agepoints = Float64[model.Tnow, model.tinit],  # Ma
        Tpoints = Float64[model.Tnow, model.Tinit],    # Degrees C
        T₀ = Float64[model.Tnow, model.Tinit],
        ΔT = Float64[model.ΔTnow, model.ΔTinit],
    )

    # Default: No unconformity is imposed
    unconf = Unconformity();

    # # Uncomment this section if you wish to impose an unconformity at any point in the record
    # # Uniform distributions from Age₀ to Age₀+ΔAge, T₀ to T₀+ΔT,
    # unconf = Unconformity(
    #     agepoints = Float64[550.0,],  # Ma
    #     Tpoints = Float64[20.0,],     # Degrees C
    #     Age₀ = Float64[500,],
    #     ΔAge = Float64[80,],
    #     T₀ = Float64[0,],
    #     ΔT = Float64[40,],
    # )
    # name *= "_unconf"

## --- Invert for maximum likelihood t-T path

    # This is where the "transdimensional" part comes in
    agepoints = zeros(Float64, model.maxpoints) # Array of fixed size to hold all optional age points
    Tpoints = zeros(Float64, model.maxpoints) # Array of fixed size to hold all optional age points

    # Fill some intermediate points to give the MCMC something to work with
    Tr = 150 # Residence temperature
    npoints = model.minpoints
    agepoints[1:npoints] .= range(0, model.tinit, npoints)
    Tpoints[1:npoints] .= Tr  # Degrees C

    # Run Markov Chain
    @time (tpointdist, Tpointdist, ndist, HeAgedist, lldist, acceptancedist, σⱼtdist, σⱼTdist) = MCMC_vartcryst(data, model, npoints, agepoints, Tpoints, boundary, unconf, detail)
    @info """tpointdist & Tpointdist collected, size: $(size(Tpointdist))
    Mean log-likelihood: $(nanmean(view(lldist, model.burnin:model.nsteps)))
    Mean acceptance rate: $(nanmean(view(acceptancedist, model.burnin:model.nsteps)))
    Mean npoints: $(nanmean(view(ndist, model.burnin:model.nsteps)))
    Mean σⱼₜ: $(nanmean(view(σⱼtdist,model.burnin:model.nsteps)))
    Mean σⱼT: $(nanmean(view(σⱼTdist, model.burnin:model.nsteps)))
    """

    # Save results using JLD
    # # Compressed:
    # using JLD
    # using JLD: @write
    # jldopen("$name.jld", "w", compress=true) do file
    #     @write file tpointdist
    #     @write file Tpointdist
    #     @write file ndist
    #     @write file HeAgedist
    #     @write file lldist
    #     @write file model
    # end
    # Uncompresed:
    # @save "$name.jld" tpointdist Tpointdist ndist HeAgedist lldist acceptancedist model
    # To read all variables from file to local workspace:
    # @load "filename.jld"

    # # Alternatively, save as MAT file
    # using MAT
    # matwrite("$name.mat", Dict(
    #     "tpointdist"=>tpointdist,
    #     "Tpointdist"=>Tpointdist,
    #     "ndist"=>ndist,
    #     "HeAgedist"=>HeAgedist,
    #     "lldist"=>lldist,
    #     "acceptancedist"=>acceptancedist,
    #     "model"=>Dict(
    #         replace.(string.(keys(model)), "σ"=>"sigma", "λ"=>"lambda", "Δ"=>"Delta") .=> values(model)
    #     )
    # ), compress=true)

    # Plot log likelihood distribution
    h = plot(lldist, xlabel="Step number", ylabel="Log likelihood", label="", framestyle=:box)
    savefig(h, name*"_lldist.pdf")
    display(h)

## --- Plot model ages vs observed ages in age-eU space (zircon)

    eU = data.U+.238*data.Th # Used only for plotting
    tz = containsi.(data.mineral, "zircon")
    if any(tz)
        h = scatter(eU[tz], data.HeAge[tz], 
            yerror=2*data.HeAge_sigma[tz], 
            label="Data (2σ)", 
            color=:black, 
            xlabel="eU (ppm)",
            ylabel="Age (Ma)",
            framestyle=:box,
        )
        zircon_agedist = HeAgedist[tz, model.burnin:end]
        m = nanmean(zircon_agedist, dims=2)
        l = nanpctile(zircon_agedist, 2.5, dims=2)
        u = nanpctile(zircon_agedist, 97.5, dims=2)
        scatter!(h, eU[tz], m, 
            yerror=(m-l, u-m), 
            label="Model + 95%CI", 
            color=:blue, 
            msc=:blue,
        )
        savefig(h, name*"_zircon_Age-eU.pdf")
        display(h)
    end

## --- Plot model ages vs observed ages in age-eU space (apatite)

    ta = containsi.(data.mineral, "apatite")
    if any(ta)
        h = scatter(eU[ta], data.HeAge[ta], 
            yerror=2*data.HeAge_sigma[ta], 
            label="Data (2σ)", 
            color=:black, 
            xlabel="eU (ppm)",
            ylabel="Age (Ma)",
            framestyle=:box,
        )
        apatite_agedist = HeAgedist[ta,model.burnin:end]
        m = nanmean(apatite_agedist, dims=2)
        l = nanpctile(apatite_agedist, 2.5, dims=2)
        u = nanpctile(apatite_agedist, 97.5, dims=2)
        scatter!(h, eU[ta], m, 
            yerror=(m-l, u-m), 
            label="Model + 95%CI", 
            color=:blue, 
            msc=:blue,
        )
        savefig(h, name*"_apatite_Age-eU.pdf")
        display(h)
    end

## --- Plot moving average of acceptance distribution

    h = plot(movmean(acceptancedist,100), label="", framestyle=:box)
    plot!(xlabel="Step number", ylabel="acceptance probability (mean of 100)")
    savefig(h, name*"_acceptance.pdf")
    display(h)

## --- Create image of paths

    # Desired rsolution of resulting image
    xresolution = 2000
    yresolution = 1000
    burnin = model.burnin

    # Resize the post-burnin part of the stationary distribution
    tTdist = Array{Float64}(undef, xresolution, model.nsteps-burnin)
    xq = range(0, model.tinit, length=xresolution)
    @time @inbounds for i = 1:model.nsteps-burnin
        linterp1s!(view(tTdist,:,i), view(tpointdist,:,i+burnin), view(Tpointdist,:,i+burnin), xq)
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

    #plot!([659, 717.4, 717.4, 659, 659],[0, 0, 650, 650, 0], fill=true, color=:white, alpha=0.6) #Sturtian glacial
    plot!([635.5, 717.4, 717.4, 635.5, 635.5],[0, 0, 650, 650, 0], fill=true, color=:white, alpha=0.5) #Sturtian & Marinoan glacial
    #plot!([635.5, 650.0, 650.0, 635.5, 635.5],[0, 0, 650, 650, 0], fill=true, color=:white, alpha=0.6) #Marinoan glacial
    #plot!([480, 640, 640, 480, 480],[0, 0, 50, 50, 0], linestyle = :dot, color=:black, linewidth=1.25) # t-T box 640 to 480 Ma, 0-50°C

    savefig(k, name*"_tT.pdf")
    display(k)

## ---

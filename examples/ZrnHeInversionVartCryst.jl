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
    using StatGeochem
    using JLD: @load, @save
    using Plots

    using Thermochron

    # Make sure we're running in the directory where the script is located
    cd(@__DIR__)

    # # # # # # # # # # Choice of regional thermochron data # # # # # # # # # #

    # Literature samples from Guenthner et al. 2013 (AJS), Minnesota
    name = "Minnesota"
    ds = importdataset("minnesota.csv", ',', importas=:Tuple);

## --- Prepare problem

    burnin = 350_000
    model = (
        nsteps = 500_000, # How many steps of the Markov chain should we run?
        burnin = burnin, # How long should we wait for MC to converge (become stationary)
        dr = 1.0,    # Radius step, in microns
        dt = 10.0,   # time step size in Myr
        dTmax = 15.0, # Maximum reheating/burial per model timestep
        TInit = 400.0, # Initial model temperature (in C) (i.e., crystallization temperature)
        ΔTInit = -50.0, # TInit can vary from TInit to TInit+ΔTinit
        TNow = 0.0, # Current surface temperature (in C)
        ΔTNow = 10.0, # TNow may vary from TNow to TNow+ΔTNow
        tInitMax = 4000.0, # Ma -- forbid anything older than this
        minPoints = 10,  # Minimum allowed number of t-T points
        maxPoints = 50, # Maximum allowed number of t-T points
        simplified = false, # Prefer simpler tT paths?
        # Diffusion parameters
        DzEa = 165.0, # kJ/mol
        DzD0 = 193188.0, # cm^2/sec
        DN17Ea = 71.0, # kJ/mol
        DN17D0 = 0.0034, #6.367E-3 # cm^2/sec
        # Model uncertainty is not well known (depends on annealing parameters,
        # decay constants, diffusion parameters, etc.), but is certainly non-zero.
        # Here we add (in quadrature) a blanket model uncertainty of 25 Ma.
        σModel = 25.0, # Ma
        σAnnealing = 35.0, # Initial annealing uncertainty [Ma]
        λAnnealing = 10 ./ burnin # Annealing decay [1/n]
    )


    # Populate data NamedTuple from imported dataset
    data = (
        halfwidth = ds.Halfwidth_um,            # Crystal half-width, in microns
        U = ds.U238_ppm,                        # U concentration, in PPM
        Th = ds.Th232_ppm,                      # Th concentration, in PPM
        HeAge = ds.HeAge_Ma_raw,                # He age, in Ma
        HeAge_sigma = ds.HeAge_Ma_sigma_raw,    # He age uncertainty (1-sigma), in Ma
        CrystAge = ds.CrystAge_Ma,              # Crystallization age, in Ma
    )

    # Sort out crystallization ages and start time
    map!(x->max(x, model.tInitMax), data.CrystAge, data.CrystAge)
    model = (model...,
        tInit = ceil(maximum(data.CrystAge)/model.dt) * model.dt,
    )

    # Boundary conditions (e.g. 10C at present and 650 C at the time of zircon formation).
    boundary = (
        agePoints = Float64[model.TNow, model.tInit],  # Ma
        TPoints = Float64[model.TNow, model.TInit],    # Degrees C
        T₀ = Float64[model.TNow, model.TInit],
        ΔT = Float64[model.ΔTNow, model.ΔTInit],
    )

    # Default: No unconformity is imposed
    unconf = (
        agePoints = Float64[],  # Ma
        TPoints = Float64[],    # Degrees C
    )

    # # Uncomment this section if you wish to impose an unconformity at any point in the record
    # # Uniform distributions from Age₀ to Age₀+ΔAge, T₀ to T₀+ΔT,
    # unconf = (;
    #     agePoints = Float64[560.0,],  # Ma
    #     TPoints = Float64[20.0,],     # Degrees C
    #     Age₀ = Float64[500,],
    #     ΔAge = Float64[80,],
    #     T₀ = Float64[0,],
    #     ΔT = Float64[50,],
    # )

    # Add additional vectors for proposed unconformity and boundary points
    unconf = (unconf...,
        agePointsₚ = similar(unconf.agePoints),
        TPointsₚ = similar(unconf.TPoints),
    )
    boundary = (boundary...,
        agePointsₚ = similar(boundary.agePoints),
        TPointsₚ = similar(boundary.TPoints),
    )
    model = (model...,
        ageSteps = Array{Float64}(model.tInit-model.dt/2 : -model.dt : 0+model.dt/2),
        tSteps = Array{Float64}(0+model.dt/2 : model.dt : model.tInit-model.dt/2),
    );

## --- Test proscribed t-T paths with Neoproterozoic exhumation step

    # # Generate T path to test
    # Tr = 150
    # T0 = 30
    # agePoints = Float64[model.tInit, model.tInit*29/30,   720, 580, 250,  0] # Age (Ma)
    # TPoints  =  Float64[model.TInit,        Tr+T0, Tr+T0,  T0,  70, 10] # Temp. (C)
    # TSteps = linterp1s(agePoints,TPoints,model.ageSteps)
    #
    # # Plot t-T path
    # plot(model.ageSteps,TSteps,xflip=true)
    #
    # # Calculate model ages
    # calcHeAges = Array{Float64}(undef, size(data.HeAge))
    # @time pr = anneal(model.dt, model.tSteps, TSteps, :zrdaam)
    # zircons = Array{Zircon{Float64}}(undef, length(data.halfwidth))
    # @time for i=1:length(zircons)
    #     # Iterate through each grain, calculate the modeled age for each
    #     first_index = 1 + floor(Int64,(model.tInit - data.CrystAge[i])/model.dt)
    #     zircons[i] = Zircon(data.halfwidth[i], model.dr, data.U[i], data.Th[i], model.dt, model.ageSteps[first_index:end])
    #     calcHeAges[i] = HeAgeSpherical(zircons[i], @views(TSteps[first_index:end]), @views(pr[first_index:end,first_index:end]), model)
    # end
    #
    # # Plot Comparison of results
    # eU = data.U+.238*data.Th # Used only for plotting
    # p2 = plot(eU, calcHeAges, seriestype=:scatter,label="Model")
    # plot!(p2, eU, data.HeAge, yerror=data.HeAge_sigma*2, seriestype=:scatter, label="Data")
    # xlabel!(p2,"eU"); ylabel!(p2,"Age (Ma)")
    # display(p2)
    #
    # # Check log likelihood
    # ll = sum(-(calcHeAges - data.HeAge).^2 ./ (2 .* data.HeAge_sigma.^2))
    # ll = sum(-(calcHeAges - data.HeAge).^2 ./ (2 .* σₐ.^2))


## --- Invert for maximum likelihood t-T path

    # This is where the "transdimensional" part comes in
    agePoints = Array{Float64}(undef, model.maxPoints) # Array of fixed size to hold all optional age points
    TPoints = Array{Float64}(undef, model.maxPoints) # Array of fixed size to hold all optional age points

    # Fill some intermediate points to give the MCMC something to work with
    Tr = 250 # Residence temperature
    nPoints = 5
    agePoints[1:nPoints] .= (model.tInit/30, model.tInit/4, model.tInit/2, model.tInit-model.tInit/4, model.tInit-model.tInit/30) # Ma
    TPoints[1:nPoints] .= Tr  # Degrees C

    # # (Optional) Initialize with something close to the expected path
    # Tr = 150
    # T0 = 30
    # nPoints = 4
    # agePoints[1:nPoints] = Float64[model.tInit*29/30,   720, 510, 250] # Age (Ma)
    # TPoints[1:nPoints]  =  Float64[            Tr+T0, Tr+T0,  T0,  70] # Temp. (C)

    # Run Markov Chain
    @time (TStepdist, HeAgedist, ndist, lldist, acceptancedist) = MCMC_vartcryst(data, model, nPoints, agePoints, TPoints, unconf, boundary)

    # # Save results using JLD
    @save string(name, ".jld") TStepdist HeAgedist lldist acceptancedist model
    # To read in from file, equivalently
    # @load "filename.jld"

    # Plot log likelihood distribution
    h = plot(lldist, xlabel="Step number", ylabel="Log likelihood", label="")
    savefig(h, name*"_lldist.pdf")
    display(h)

## ---  Plot sample age-eU correlations

    eU = data.U+.238*data.Th # Used only for plotting
    h = scatter(eU,HeAgedist[:,model.burnin:model.burnin+50], label="")
    plot!(h, eU, data.HeAge, yerror=data.HeAge_sigma, seriestype=:scatter, label="Data")
    xlabel!(h,"eU (ppm)"); ylabel!(h,"Age (Ma)")
    savefig(h, name*"_Age-eU.pdf")
    display(h)


## --- Create image of paths

    # Desired rsolution of resulting image
    xresolution = 2000
    yresolution = 1000

    # Resize the post-burnin part of the stationary distribution
    tTdist = Array{Float64}(undef, xresolution, size(TStepdist,2)-model.burnin)
    xq = range(0,model.tInit,length=xresolution)
    @time @inbounds for i = 1:size(TStepdist,2)-model.burnin
        linterp1!(view(tTdist,:,i), model.tSteps, view(TStepdist,:,i+model.burnin), xq)
    end

    # Calculate composite image
    ybinedges = range(model.TNow, model.TInit, length=yresolution+1)
    tTimage = zeros(yresolution, size(tTdist,1))
    @time @inbounds for i=1:size(tTdist,1)
        histcounts!(view(tTimage,:,i), view(tTdist,i,:), ybinedges)
    end

## --- Plot image with 'ylcn' custom colorscale

    # Prepare axes
    k = plot(layout = grid(1,2, widths=[0.94, 0.06]), framestyle=:box)

    # Plot image with colorscale in first subplot
    A = imsc(reverse(tTimage,dims=2), ylcn, 0, nanpctile(tTimage[:],98.5))
    plot!(k[1], xlabel="Time (Ma)",ylabel="Temperature (°C)",yticks=0:50:400,xticks=0:500:3500,yminorticks=5,xminorticks=5,tick_dir=:out,framestyle=:box)
    plot!(k[1],xq,cntr(ybinedges),A,yflip=true,xflip=true,legend=false,aspectratio=3500/400/1.5,xlims=(0,3500),ylims=(0,400))

    # Add colorbar in second subplot
    cb = imsc(repeat(0:100, 1, 10), ylcn, 0, 100)
    plot!(k[2], 0:0.01:0.1, 0:0.01:1, cb, ylims=(0,1), xticks=false, framestyle=:box, yflip=false)

    #plot!([659, 717.4, 717.4, 659, 659],[0, 0, 650, 650, 0], fill=true, color=:white, alpha=0.6) #Sturtian glacial
    plot!([635.5, 717.4, 717.4, 635.5, 635.5],[0, 0, 650, 650, 0], fill=true, color=:white, alpha=0.5) #Sturtian & Marinoan glacial
    #plot!([635.5, 650.0, 650.0, 635.5, 635.5],[0, 0, 650, 650, 0], fill=true, color=:white, alpha=0.6) #Marinoan glacial
    #plot!([480, 640, 640, 480, 480],[0, 0, 50, 50, 0], linestyle = :dot, color=:black, linewidth=1.25) # t-T box 640 to 480 Ma, 0-50°C

    savefig(k, name*"mcmc-zrdaam.pdf")
    display(k)

## ---

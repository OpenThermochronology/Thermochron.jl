# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                            ZrnHeInversion.jl                                  #
#                                                                               #
#       A transdimensional Bayesian Markov-chain Monte Carlo code for           #
#   inverting time-temperature paths from zircon U-Th/He ages, using a 1D       #
#   Crank-Nicholson finite difference diffusion model along with the zircon     #
#   damage and annealing model "ZRDAAM" from Guenthner et al, 2013 (AJS).       #
#   Requires uncorrected ZHe ages, grain size, [U] and [Th] in ppm.             #
#                                                                               #
#       This file contains the MCMC code, while "ZrnHe.jl" contains the damage  #
#   annealing function and the Crank-Nicholson code.                            #
#                                                                               #
#   © 2022 C. Brenhin Keller and Kalin McDannell                                #                                                                            #
#                                                                               #
#   If running for the first time, make sure you have instantiated the          #
#   manifest that came with this file                                           #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
## ---  Load required packages
    using Statistics, LinearAlgebra
    using StatsBase: fit, Histogram
    using ProgressMeter: @showprogress
    using StatGeochem
    using JLD: @load, @save
    using Plots; gr() # Use the GR backend for plotting

    using Thermochron

    # Make sure we're running in the directory where the script is located
    cd(@__DIR__)

    # # # # # # # # # # Choice of regional thermochron data # # # # # # # # # #

    # Literature samples from Guenthner et al. 2013 (AJS), Minnesota
    name = "MinnesotaInversion"
    ds = importdataset("minnesota.csv", ',', importas=:Tuple)


## --- Prepare problem

    burnin = 100000
    model = (
        nsteps = 200000, # How many steps of the Markov chain should we run?
        burnin = burnin, # How long should we wait for MC to converge (become stationary)
        dr = 1.0,    # Radius step, in microns
        dt = 10.0,   # time step size in Myr
        dTmax = 10.0, # Maximum reheating/burial per model timestep (to prevent numerical overflow)
        TInit = 400.0, # Initial model temperature (in C) (i.e., crystallization temperature)
        ΔTInit = -50.0, # TInit can vary from TInit to TInit+ΔTinit
        TNow = 0.0, # Current surface temperature (in C)
        ΔTNow = 10.0, # TNow may vary from TNow to TNow+ΔTNow
        tInitMax = 4000.0, # Ma -- forbid anything older than this
        maxPoints = 40, # Maximum allowed number of t-T points
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

    tSteps = Array{Float64,1}(0+model.dt/2 : model.dt : model.tInit-model.dt/2)
    ageSteps = reverse(tSteps)

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

## --- Test proscribed t-T paths with Neoproterozoic exhumation step

    # # Generate T path to test
    # Tr = 150
    # T0 = 30
    # agePoints = Float64[model.tInit, model.tInit*29/30,   720, 580, 250,  0] # Age (Ma)
    # TPoints  =  Float64[model.TInit,        Tr+T0, Tr+T0,  T0,  70, 10] # Temp. (C)
    # TSteps = linterp1s(agePoints,TPoints,ageSteps)
    #
    # # Plot t-T path
    # plot(ageSteps,TSteps,xflip=true)
    #
    # # Calculate model ages
    # calcHeAges = Array{Float64}(undef, size(data.HeAge))
    # @time pr = anneal(model.dt, tSteps, TSteps, :zrdaam)
    # zircons = Array{Zircon{Float64}}(undef, length(data.halfwidth))
    # @time for i=1:length(zircons)
    #     # Iterate through each grain, calculate the modeled age for each
    #     first_index = 1 + floor(Int64,(model.tInit - data.CrystAge[i])/model.dt)
    #     zircons[i] = Zircon(data.halfwidth[i], model.dr, data.U[i], data.Th[i], model.dt, ageSteps[first_index:end])
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

    # Add additional vectors for proposed unconformity and boundary points
    unconf = (unconf...,
        agePointsₚ = similar(unconf.agePoints),
        TPointsₚ = similar(unconf.TPoints),
    )
    boundary = (boundary...,
        agePointsₚ = similar(boundary.agePoints),
        TPointsₚ = similar(boundary.TPoints),
    )

    # This is where the "transdimensional" part comes in
    agePoints = Array{Float64}(undef, model.maxPoints+1) # Array of fixed size to hold all optional age points
    TPoints = Array{Float64}(undef, model.maxPoints+1) # Array of fixed size to hold all optional age points
    # Fill some intermediate points to give the MCMC something to work with
    Tr = 250 # Residence temperature
    nPoints = 5
    agePoints[1:nPoints] .= [model.tInit/30,model.tInit/4,model.tInit/2,model.tInit-model.tInit/4,model.tInit-model.tInit/30] # Ma
    TPoints[1:nPoints] .+ Tr  # Degrees C

    # # (Optional) Initialize with something close to the expected path
    # Tr = 150
    # T0 = 30
    # nPoints = 4
    # agePoints[1:nPoints] = Float64[model.tInit*29/30,   720, 510, 250]) # Age (Ma)
    # TPoints[1:nPoints]  =  Float64[       Tr+T0, Tr+T0,  T0,  70]) # Temp. (C)

    function MCMC_vartcryst(data, model, nPoints, agePoints, TPoints, unconf, boundary)
        # Calculate model ages for initial proposal
        TSteps = linterp1s([view(agePoints, 1:nPoints) ; boundary.agePoints ; unconf.agePoints],
                           [view(TPoints, 1:nPoints) ; boundary.TPoints ; unconf.TPoints], ageSteps)
        calcHeAges = Array{Float64}(undef, size(data.HeAge))
        pr = anneal(model.dt, tSteps, TSteps, :zrdaam) # Damage annealing history

        # Prepare a Mineral object for each analysis
        zircons = Array{Zircon{Float64}}(undef, length(data.halfwidth))
        for i=1:length(zircons)
            # Iterate through each grain, calculate the modeled age for each
            first_index = 1 + floor(Int64,(model.tInit - data.CrystAge[i])/model.dt)
            zircons[i] = Zircon(data.halfwidth[i], model.dr, data.U[i], data.Th[i], model.dt, ageSteps[first_index:end])
            calcHeAges[i] = HeAgeSpherical(zircons[i], @views(TSteps[first_index:end]), @views(pr[first_index:end,first_index:end]), model)
        end

        # Simulated annealing of uncertainty
        σₐ = simannealsigma.(1, data.HeAge_sigma; params=model)
        σₙ = simannealsigma.(model.nsteps, data.HeAge_sigma; params=model)

        # Log-likelihood for initial proposal
        ll = normpdf_ll(data.HeAge, σₐ, calcHeAges)
        if model.simplified
            ll -= log(nPoints)
        end

        # Variables to hold proposals
        llₚ = ll
        nPointsₚ = nPoints
        agePointsₚ = similar(agePoints)
        TPointsₚ = similar(TPoints)
        TStepsₚ = similar(TSteps)
        calcHeAgesₚ = similar(calcHeAges)

        # Distributions to populate
        HeAgeDist = Array{Float64}(undef, length(data.HeAge), model.nsteps)
        TStepDist = Array{Float64}(undef, length(tSteps), model.nsteps)
        llDist = Array{Float64}(undef, model.nsteps)
        nDist = zeros(Int, model.nsteps)
        acceptanceDist = zeros(Bool, model.nsteps)

        # Standard deviations of Gaussian proposal distributions for temperature and time
        t_sigma = model.tInit/60
        T_sigma = model.TInit/60

        # Proposal probabilities (must sum to 1)
        move = 0.65
        birth = 0.15
        death = 0.15 # Should equal birth
        boundarymove = 0.05
        maxattempts = 1000

        @showprogress 10 "Running MCMC..." for n=1:model.nsteps

            # Copy proposal from last accepted solution
            nPointsₚ = nPoints
            copyto!(agePointsₚ, agePoints)
            copyto!(TPointsₚ, TPoints)
            copyto!(TStepsₚ, TSteps)
            copyto!(unconf.agePointsₚ, unconf.agePoints)
            copyto!(unconf.TPointsₚ, unconf.TPoints)
            copyto!(boundary.TPointsₚ, boundary.TPoints)

            # Adjust the proposal
            r = rand()
            if r < move
                # Move the age of one model point
                for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                    k = ceil(Int, rand() * nPoints)

                    agePointsₚ[k] += randn() * t_sigma
                    if agePointsₚ[k] < model.dt
                        # Don't let any point get too close to 0
                        agePointsₚ[k] += (model.dt - agePointsₚ[k])
                    elseif agePointsₚ[k] > (model.tInit - model.dt)
                        # Don't let any point get too close to tInit
                        agePointsₚ[k] -= (agePointsₚ[k] - (model.tInit - model.dt))
                    end
                    # Move the Temperature of one model point
                    if TPointsₚ[k] < 0
                        # Don't allow T<0
                        TPointsₚ[k] = 0
                    elseif TPointsₚ[k] > model.TInit
                        # Don't allow T>TInit
                        TPointsₚ[k] = model.TInit
                    end

                    # Interpolate proposed t-T path
                    TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                        [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], ageSteps)

                    # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                    maximum(diff(TStepsₚ)) < model.dTmax && break

                    # Copy last accepted solution to re-modify if we don't break
                    copyto!(agePointsₚ, agePoints)
                    copyto!(TPointsₚ, TPoints)
                end
            elseif (r < move+birth) && (nPointsₚ < model.maxPoints)
                # Birth: add a new model point
                nPointsₚ += 1
                for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                    agePointsₚ[nPointsₚ] = rand()*model.tInit
                    TPointsₚ[nPointsₚ] = rand()*model.TInit

                    # Interpolate proposed t-T path
                    TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                        [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], ageSteps)

                    # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                    maximum(diff(TStepsₚ)) < model.dTmax && break
                end
            elseif (r < move+birth+death) && (r >= move+birth) && (nPointsₚ > 1)
                # Death: remove a model point
                nPointsₚ -= 1 # Delete last point in array from proposal
                for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                    k = ceil(Int, rand()*nPoints) # Choose point to delete
                    agePointsₚ[k] = agePointsₚ[nPoints]
                    TPointsₚ[k] = TPointsₚ[nPoints]

                    # Interpolate proposed t-T path
                    TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                        [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], ageSteps)

                    # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                    maximum(diff(TStepsₚ)) < model.dTmax && break
                end
            else
                # Move boundary conditions
                for i=1:maxattempts # Try maxattempts times to satisfy the reheating rate limit
                    # Move the temperatures of the starting and ending boundaries
                    @. boundary.TPointsₚ = boundary.T₀ + rand()*boundary.ΔT

                    # If there's an imposed unconformity, adjust within parameters
                    if length(unconf.agePoints) > 0
                        @. unconf.agePointsₚ = unconf.Age₀ + rand()*unconf.ΔAge
                        @. unconf.TPointsₚ = unconf.T₀ + rand()*unconf.ΔT
                    end


                    # Recalculate interpolated proposed t-T path
                    TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                            [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], ageSteps)

                    # Accept the proposal (and break out of the loop) if it satisfies the maximum reheating rate
                    maximum(diff(TStepsₚ)) < model.dTmax && break

                    # Copy last accepted solution to re-modify if we don't break
                    copyto!(unconf.agePointsₚ, unconf.agePoints)
                    copyto!(unconf.TPointsₚ, unconf.TPoints)
                    copyto!(boundary.TPointsₚ, boundary.TPoints)
                end
            end

            # Recalculate interpolated proposed t-T path
            TStepsₚ = linterp1s([view(agePointsₚ, 1:nPointsₚ) ; boundary.agePoints ; unconf.agePointsₚ],
                                [view(TPointsₚ, 1:nPointsₚ) ; boundary.TPointsₚ ; unconf.TPointsₚ], ageSteps)

             # Calculate model ages for each grain
            anneal!(pr, model.dt, tSteps, TStepsₚ, :zrdaam)
            for i=1:length(zircons)
                first_index = 1 + floor(Int64,(model.tInit - data.CrystAge[i])/model.dt)
                calcHeAgesₚ[i] = HeAgeSpherical(zircons[i], @views(TStepsₚ[first_index:end]), @views(pr[first_index:end,first_index:end]), model)
            end

            # Calculate log likelihood of proposal
            σₐ .= simannealsigma.(n, data.HeAge_sigma; params=model)
            llₚ = normpdf_ll(data.HeAge, σₐ, calcHeAgesₚ)
            llₗ = normpdf_ll(data.HeAge, σₐ, calcHeAges) # Recalulate last one too with new σₐ
            if model.simplified # slightly penalize more complex t-T paths
                llₚ -= log(nPointsₚ)
                llₗ -= log(nPoints)
            end

            # Accept or reject proposal based on likelihood
            # To avoid numerical problems with diffusion code, also reject proposal
            # if maximum proposed heating rate is greater than 25C per timestep.
            # (Fast cooling should not be a problem, however)
            if log(rand()) < (llₚ - llₗ)
                ll = llₚ
                nPoints = nPointsₚ
                copyto!(agePoints, agePointsₚ)
                copyto!(TPoints, TPointsₚ)
                copyto!(unconf.agePoints, unconf.agePointsₚ)
                copyto!(unconf.TPoints, unconf.TPointsₚ)
                copyto!(boundary.TPoints, boundary.TPointsₚ)
                copyto!(calcHeAges, calcHeAgesₚ)

                # These are saved for ouput, but not critical to the function of the MCMC loop
                copyto!(TSteps, TStepsₚ)
                acceptanceDist[n] = true
            end

            # Record results for analysis and troubleshooting
            llDist[n] = normpdf_ll(data.HeAge, σₙ, calcHeAges) # Recalculated to constant baseline
            nDist[n] = nPoints # Distribution of # of points
            HeAgeDist[:,n] = calcHeAges # Distribution of He ages

            # This is the actual output we want -- the distribution of t-T paths (t path is always identical)
            TStepDist[:,n] = TSteps # Distribution of T paths

        end
        return (TStepDist, HeAgeDist, nDist, llDist, acceptanceDist)
    end

    # Run Markov Chain
    @time (TStepDist, HeAgeDist, nDist, llDist, acceptanceDist) = MCMC_vartcryst(data, model, nPoints, agePoints, TPoints, unconf, boundary)

    # # Save results using JLD
    @save string(name, ".jld") ageSteps tSteps TStepDist model

## ---  Plot sample age-eU correlations

    eU = data.U+.238*data.Th # Used only for plotting
    h = scatter(eU,HeAgeDist[:,model.burnin:model.burnin+50], label="")
    plot!(h, eU, data.HeAge, yerror=data.HeAge_sigma, seriestype=:scatter, label="Data")
    xlabel!(h,"eU (ppm)"); ylabel!(h,"Age (Ma)")
    savefig(h,string(name,"Age-eU.pdf"))


## --- Create image of paths

    # Resize the post-burnin part of the stationary distribution
    TStepDistResized = Array{Float64}(undef, 2001, size(TStepDist,2)-model.burnin)
    xq = collect(range(0,model.tInit,length=2001))
    for i=1:size(TStepDist,2)-model.burnin
        TStepDistResized[:,i] = linterp1s(tSteps,TStepDist[:,i+model.burnin],xq)
    end

    # Calculate composite image
    tTimage = zeros(ceil(Int, TInit)*2, size(TStepDistResized,1))
    yq = collect(0:0.5:TInit)
    for i=1:size(TStepDistResized,1)
        hist = fit(Histogram,TStepDistResized[i,:],yq,closed=:right)
        tTimage[:,i] = hist.weights
    end

## --- Plot image with 'ylcn' custom colorscale
    # Prepare axes
    k = plot(layout = grid(1,2, widths=[0.94, 0.06]), framestyle=:box)

    # Plot image with colorscale in first subplot
    A = imsc(reverse(tTimage,dims=2), ylcn, 0, nanpctile(tTimage[:],98.5))
    plot!(k[1], xlabel="Time (Ma)",ylabel="Temperature (°C)",yticks=0:50:400,xticks=0:500:3500,yminorticks=5,xminorticks=5,tick_dir=:out,framestyle=:box)
    plot!(k[1],xq,cntr(yq),A,yflip=true,xflip=true,legend=false,aspectratio=3500/400/1.5,xlims=(0,3500),ylims=(0,400))

    # Add colorbar in second subplot
    cb = imsc(repeat(0:100, 1, 10), ylcn, 0, 100)
    plot!(k[2], 0:0.01:0.1, 0:0.01:1, cb, xticks=false, framestyle=:box, yflip=false)

    #plot!([659, 717.4, 717.4, 659, 659],[0, 0, 650, 650, 0], fill=true, color=:white, alpha=0.6) #Sturtian glacial
    plot!([635.5, 717.4, 717.4, 635.5, 635.5],[0, 0, 650, 650, 0], fill=true, color=:white, alpha=0.5) #Sturtian & Marinoan glacial
    #plot!([635.5, 650.0, 650.0, 635.5, 635.5],[0, 0, 650, 650, 0], fill=true, color=:white, alpha=0.6) #Marinoan glacial
    #plot!([480, 640, 640, 480, 480],[0, 0, 50, 50, 0], linestyle = :dot, color=:black, linewidth=1.25) # t-T box 640 to 480 Ma, 0-50°C


    savefig(k,"minnesota-mcmc-zrdaam-unc.pdf")
    display(k)

    #img = plot(reverse(xq), cntr(yq), imsc(reverse(tTimage,dims=2), viridis, 0, nanpctile(tTimage[:],97.5)), yflip=true, xflip=true, legend=false,
    #     xlabel="Age (Ma)", ylabel="Temperature (°C)")
    #savefig(img,string(name,"minnesota-mcmc.pdf"))
    #display(img)

## ---
